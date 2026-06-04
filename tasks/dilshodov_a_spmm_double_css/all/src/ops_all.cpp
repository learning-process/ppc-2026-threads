#include "dilshodov_a_spmm_double_css/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "dilshodov_a_spmm_double_css/common/include/common.hpp"

namespace dilshodov_a_spmm_double_css {

namespace {
constexpr double kEps = 1e-10;

bool HasValidDimensions(const SparseMatrixCCS &m) {
  return m.rows_count > 0 && m.cols_count > 0;
}

bool HasValidContainers(const SparseMatrixCCS &m) {
  if (m.col_ptrs.size() != static_cast<size_t>(m.cols_count) + 1) {
    return false;
  }
  if (m.row_indices.size() != m.values.size()) {
    return false;
  }
  if (m.col_ptrs.empty() || m.col_ptrs.front() != 0) {
    return false;
  }
  if (m.col_ptrs.back() < 0) {
    return false;
  }
  if (static_cast<size_t>(m.col_ptrs.back()) != m.values.size()) {
    return false;
  }
  return true;
}

bool HasValidColumnOrdering(const SparseMatrixCCS &m) {
  for (int j = 0; j < m.cols_count; ++j) {
    if (m.col_ptrs[j] > m.col_ptrs[j + 1]) {
      return false;
    }
    int prev_row = -1;
    for (int idx = m.col_ptrs[j]; idx < m.col_ptrs[j + 1]; ++idx) {
      const int row = m.row_indices[idx];
      if (row < 0 || row >= m.rows_count) {
        return false;
      }
      if (row <= prev_row) {
        return false;
      }
      prev_row = row;
    }
  }
  return true;
}

bool IsValidCCS(const SparseMatrixCCS &m) {
  return HasValidDimensions(m) && HasValidContainers(m) && HasValidColumnOrdering(m);
}

void AccumulateColumnProduct(const SparseMatrixCCS &lhs, const SparseMatrixCCS &rhs, int rhs_col,
                             std::vector<double> &acc, std::vector<int> &marker, std::vector<int> &used_rows,
                             std::vector<std::pair<int, double>> &col) {
  used_rows.clear();
  col.clear();

  for (int idx_b = rhs.col_ptrs[rhs_col]; idx_b < rhs.col_ptrs[rhs_col + 1]; ++idx_b) {
    const int k = rhs.row_indices[idx_b];
    const double rhs_value = rhs.values[idx_b];

    for (int idx_a = lhs.col_ptrs[k]; idx_a < lhs.col_ptrs[k + 1]; ++idx_a) {
      const int row = lhs.row_indices[idx_a];
      const double value = lhs.values[idx_a] * rhs_value;

      if (marker[row] != rhs_col) {
        marker[row] = rhs_col;
        acc[row] = value;
        used_rows.push_back(row);
      } else {
        acc[row] += value;
      }
    }
  }

  col.reserve(used_rows.size());
  for (int row : used_rows) {
    const double value = acc[row];
    if (std::abs(value) > kEps) {
      col.emplace_back(row, value);
    }
  }

  std::ranges::sort(col, {}, &std::pair<int, double>::first);
}

void ComputeLocalSlab(const SparseMatrixCCS &lhs, const SparseMatrixCCS &rhs, int j_start, int j_end,
                      std::vector<int> &local_col_ptrs, std::vector<int> &local_row_indices,
                      std::vector<double> &local_values) {
  const int local_cols = j_end - j_start;
  std::vector<std::vector<std::pair<int, double>>> column_results(local_cols);

#pragma omp parallel default(none) shared(lhs, rhs, column_results, j_start, j_end)
  {
    std::vector<double> acc(lhs.rows_count, 0.0);
    std::vector<int> marker(lhs.rows_count, -1);
    std::vector<int> used_rows;
    used_rows.reserve(256);

#pragma omp for schedule(dynamic)
    for (int j = j_start; j < j_end; ++j) {
      AccumulateColumnProduct(lhs, rhs, j, acc, marker, used_rows, column_results[j - j_start]);
    }
  }

  local_col_ptrs.assign(static_cast<size_t>(local_cols) + 1, 0);
  for (int col_idx = 0; col_idx < local_cols; ++col_idx) {
    local_col_ptrs[col_idx + 1] = local_col_ptrs[col_idx] + static_cast<int>(column_results[col_idx].size());
  }
  const int local_nnz = local_col_ptrs[local_cols];
  local_row_indices.resize(local_nnz);
  local_values.resize(local_nnz);

  for (int col_idx = 0; col_idx < local_cols; ++col_idx) {
    int write = local_col_ptrs[col_idx];
    for (const auto &[row, value] : column_results[col_idx]) {
      local_row_indices[write] = row;
      local_values[write] = value;
      ++write;
    }
  }
}

}  // namespace

DilshodovASpmmDoubleCssAll::DilshodovASpmmDoubleCssAll(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool DilshodovASpmmDoubleCssAll::ValidationImpl() {
  const auto &lhs = std::get<0>(GetInput());
  const auto &rhs = std::get<1>(GetInput());
  return IsValidCCS(lhs) && IsValidCCS(rhs) && lhs.cols_count == rhs.rows_count;
}

bool DilshodovASpmmDoubleCssAll::PreProcessingImpl() {
  GetOutput() = SparseMatrixCCS{};
  return true;
}

bool DilshodovASpmmDoubleCssAll::RunImpl() {
  const auto &lhs = std::get<0>(GetInput());
  const auto &rhs = std::get<1>(GetInput());
  auto &out = GetOutput();

  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  out.rows_count = lhs.rows_count;
  out.cols_count = rhs.cols_count;
  out.col_ptrs.assign(static_cast<size_t>(out.cols_count) + 1, 0);
  out.row_indices.clear();
  out.values.clear();

  std::vector<int> col_starts(static_cast<size_t>(size) + 1);
  for (int proc = 0; proc <= size; ++proc) {
    col_starts[proc] = static_cast<int>((static_cast<int64_t>(proc) * rhs.cols_count) / size);
  }

  const int j_start = col_starts[rank];
  const int j_end = col_starts[rank + 1];
  const int local_cols = j_end - j_start;

  std::vector<int> local_col_ptrs;
  std::vector<int> local_row_indices;
  std::vector<double> local_values;

  if (local_cols > 0) {
    ComputeLocalSlab(lhs, rhs, j_start, j_end, local_col_ptrs, local_row_indices, local_values);
  } else {
    local_col_ptrs.assign(1, 0);
  }

  const int local_nnz = local_col_ptrs.empty() ? 0 : local_col_ptrs.back();

  std::vector<int> nnz_per_proc(size, 0);
  MPI_Allgather(&local_nnz, 1, MPI_INT, nnz_per_proc.data(), 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<int> nnz_displs(size, 0);
  int total_nnz = 0;
  for (int proc = 0; proc < size; ++proc) {
    nnz_displs[proc] = total_nnz;
    total_nnz += nnz_per_proc[proc];
  }

  out.row_indices.resize(total_nnz);
  out.values.resize(total_nnz);

  MPI_Allgatherv(local_row_indices.data(), local_nnz, MPI_INT, out.row_indices.data(), nnz_per_proc.data(),
                 nnz_displs.data(), MPI_INT, MPI_COMM_WORLD);
  MPI_Allgatherv(local_values.data(), local_nnz, MPI_DOUBLE, out.values.data(), nnz_per_proc.data(), nnz_displs.data(),
                 MPI_DOUBLE, MPI_COMM_WORLD);

  std::vector<int> cols_per_proc(size, 0);
  std::vector<int> col_displs(size, 0);
  for (int proc = 0; proc < size; ++proc) {
    cols_per_proc[proc] = col_starts[proc + 1] - col_starts[proc];
    col_displs[proc] = col_starts[proc];
  }

  std::vector<int> local_inner(static_cast<size_t>(local_cols));
  for (int col_idx = 0; col_idx < local_cols; ++col_idx) {
    local_inner[col_idx] = local_col_ptrs[col_idx] + nnz_displs[rank];
  }

  MPI_Allgatherv(local_inner.data(), local_cols, MPI_INT, out.col_ptrs.data(), cols_per_proc.data(), col_displs.data(),
                 MPI_INT, MPI_COMM_WORLD);
  out.col_ptrs[out.cols_count] = total_nnz;

  out.non_zeros = total_nnz;
  return true;
}

bool DilshodovASpmmDoubleCssAll::PostProcessingImpl() {
  GetOutput().non_zeros = static_cast<int>(GetOutput().values.size());
  return true;
}

}  // namespace dilshodov_a_spmm_double_css
