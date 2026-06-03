#include "kotelnikova_a_double_matr_mult/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "kotelnikova_a_double_matr_mult/common/include/common.hpp"

namespace kotelnikova_a_double_matr_mult {

KotelnikovaATaskALL::KotelnikovaATaskALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = SparseMatrixCCS();
}

bool KotelnikovaATaskALL::IsMatrixValid(const SparseMatrixCCS &matrix) {
  if (matrix.rows < 0 || matrix.cols < 0) {
    return false;
  }
  if (matrix.col_ptrs.size() != static_cast<size_t>(matrix.cols) + 1) {
    return false;
  }
  if (matrix.values.size() != matrix.row_indices.size()) {
    return false;
  }

  if (matrix.col_ptrs.empty() || matrix.col_ptrs[0] != 0) {
    return false;
  }

  const int total_elements = static_cast<int>(matrix.values.size());
  if (matrix.col_ptrs[matrix.cols] != total_elements) {
    return false;
  }

  for (size_t i = 0; i < matrix.col_ptrs.size() - 1; ++i) {
    if (matrix.col_ptrs[i] > matrix.col_ptrs[i + 1] || matrix.col_ptrs[i] < 0) {
      return false;
    }
  }

  for (size_t i = 0; i < matrix.row_indices.size(); ++i) {
    if (matrix.row_indices[i] < 0 || matrix.row_indices[i] >= matrix.rows) {
      return false;
    }
  }

  return true;
}

bool KotelnikovaATaskALL::ValidationImpl() {
  const auto &[a, b] = GetInput();

  if (!IsMatrixValid(a) || !IsMatrixValid(b)) {
    return false;
  }
  if (a.cols != b.rows) {
    return false;
  }

  return true;
}

bool KotelnikovaATaskALL::PreProcessingImpl() {
  const auto &[a, b] = GetInput();
  GetOutput() = SparseMatrixCCS(a.rows, b.cols);
  return true;
}

namespace {
void ComputeColumnBlock(const SparseMatrixCCS &a, const SparseMatrixCCS &b, int start_col, int end_col,
                        std::vector<std::vector<double>> &temp_columns, std::vector<int> &nnz_counts) {
  const double epsilon = 1e-10;

#pragma omp parallel for default(none) shared(a, b, start_col, end_col, temp_columns, nnz_counts, epsilon) \
    schedule(dynamic, 4)
  for (int j = start_col; j < end_col; ++j) {
    std::vector<double> &temp = temp_columns[j - start_col];
    temp.assign(a.rows, 0.0);

    for (int b_idx = b.col_ptrs[j]; b_idx < b.col_ptrs[j + 1]; ++b_idx) {
      const int k = b.row_indices[b_idx];
      const double b_val = b.values[b_idx];

      for (int a_idx = a.col_ptrs[k]; a_idx < a.col_ptrs[k + 1]; ++a_idx) {
        const int i = a.row_indices[a_idx];
        temp[i] += a.values[a_idx] * b_val;
      }
    }

    int count = 0;
    for (double val : temp) {
      if (std::abs(val) > epsilon) {
        ++count;
      }
    }
    nnz_counts[j - start_col] = count;
  }
}

void BuildLocalResult(const SparseMatrixCCS &a, int start_col, int end_col,
                      const std::vector<std::vector<double>> &temp_columns, const std::vector<int> &nnz_counts,
                      std::vector<double> &local_values, std::vector<int> &local_row_indices,
                      std::vector<int> &local_col_ptrs) {
  (void)a;

  const double epsilon = 1e-10;
  const int local_cols = end_col - start_col;

  local_col_ptrs.resize(local_cols + 1, 0);
  for (int j = 0; j < local_cols; ++j) {
    local_col_ptrs[j + 1] = local_col_ptrs[j] + nnz_counts[j];
  }

  const int total_nnz = local_col_ptrs[local_cols];
  local_values.resize(total_nnz);
  local_row_indices.resize(total_nnz);

#pragma omp parallel for default(none) shared(a, start_col, end_col, temp_columns, local_values, local_row_indices, \
                                                  local_col_ptrs, epsilon, local_cols) schedule(dynamic, 4)
  for (int j = 0; j < local_cols; ++j) {
    const std::vector<double> &temp = temp_columns[j];
    int pos = local_col_ptrs[j];
    for (size_t i = 0; i < temp.size(); ++i) {
      if (std::abs(temp[i]) > epsilon) {
        local_row_indices[pos] = static_cast<int>(i);
        local_values[pos] = temp[i];
        ++pos;
      }
    }
  }
}

void GatherAndBroadcastResult(int rank, int size, const std::vector<double> &local_values,
                              const std::vector<int> &local_row_indices, const std::vector<int> &local_col_ptrs,
                              int total_cols, int rows, SparseMatrixCCS &result) {
  std::vector<int> recv_counts(size, 0);
  std::vector<int> recv_offsets(size, 0);

  int local_nnz = static_cast<int>(local_values.size());

  MPI_Gather(&local_nnz, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<double> global_values;
  std::vector<int> global_row_indices;
  std::vector<int> global_col_ptrs;
  std::vector<int> global_col_counts(total_cols, 0);

  if (rank == 0) {
    int total_nnz = 0;
    for (int i = 0; i < size; ++i) {
      recv_offsets[i] = total_nnz;
      total_nnz += recv_counts[i];
    }
    global_values.resize(total_nnz);
    global_row_indices.resize(total_nnz);
    global_col_ptrs.resize(total_cols + 1, 0);
  }

  MPI_Gatherv(local_values.data(), local_nnz, MPI_DOUBLE, global_values.data(), recv_counts.data(), recv_offsets.data(),
              MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Gatherv(local_row_indices.data(), local_nnz, MPI_INT, global_row_indices.data(), recv_counts.data(),
              recv_offsets.data(), MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<int> local_col_counts(total_cols, 0);
  int cols_per_proc = total_cols / size;
  int remainder = total_cols % size;

  int start = (rank * cols_per_proc) + std::min(rank, remainder);
  int end = start + cols_per_proc + (rank < remainder ? 1 : 0);

  for (int j = start; j < end; ++j) {
    local_col_counts[j] = (j - start + 1 < static_cast<int>(local_col_ptrs.size()))
                              ? (local_col_ptrs[j - start + 1] - local_col_ptrs[j - start])
                              : 0;
  }

  std::vector<int> global_col_counts_tmp(total_cols, 0);
  MPI_Reduce(local_col_counts.data(), global_col_counts_tmp.data(), total_cols, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    global_col_ptrs[0] = 0;
    for (int j = 0; j < total_cols; ++j) {
      global_col_ptrs[j + 1] = global_col_ptrs[j] + global_col_counts_tmp[j];
    }

    result.values = global_values;
    result.row_indices = global_row_indices;
    result.col_ptrs = global_col_ptrs;
    result.rows = rows;
    result.cols = total_cols;
  }

  int result_values_size = 0;
  int result_row_indices_size = 0;
  int result_col_ptrs_size = 0;

  if (rank == 0) {
    result_values_size = static_cast<int>(result.values.size());
    result_row_indices_size = static_cast<int>(result.row_indices.size());
    result_col_ptrs_size = static_cast<int>(result.col_ptrs.size());
  }

  MPI_Bcast(&result_values_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&result_row_indices_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&result_col_ptrs_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    result.values.resize(result_values_size);
    result.row_indices.resize(result_row_indices_size);
    result.col_ptrs.resize(result_col_ptrs_size);
    result.rows = rows;
    result.cols = total_cols;
  }

  MPI_Bcast(result.values.data(), result_values_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(result.row_indices.data(), result_row_indices_size, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(result.col_ptrs.data(), result_col_ptrs_size, MPI_INT, 0, MPI_COMM_WORLD);
}

}  // namespace

SparseMatrixCCS KotelnikovaATaskALL::MultiplyMatricesMPIOMP(const SparseMatrixCCS &a, const SparseMatrixCCS &b) {
  int rank = -1;
  int size = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const int total_cols = b.cols;
  const int cols_per_proc = total_cols / size;
  const int remainder = total_cols % size;

  int start_col = (rank * cols_per_proc) + std::min(rank, remainder);
  int end_col = start_col + cols_per_proc + (rank < remainder ? 1 : 0);

  const int local_cols = end_col - start_col;

  if (local_cols > 0) {
    std::vector<std::vector<double>> temp_columns(local_cols);
    std::vector<int> nnz_counts(local_cols, 0);

    ComputeColumnBlock(a, b, start_col, end_col, temp_columns, nnz_counts);

    std::vector<double> local_values;
    std::vector<int> local_row_indices;
    std::vector<int> local_col_ptrs;

    BuildLocalResult(a, start_col, end_col, temp_columns, nnz_counts, local_values, local_row_indices, local_col_ptrs);

    SparseMatrixCCS result(a.rows, total_cols);
    GatherAndBroadcastResult(rank, size, local_values, local_row_indices, local_col_ptrs, total_cols, a.rows, result);
    return result;
  }
  SparseMatrixCCS result(a.rows, total_cols);
  std::vector<double> empty_values;
  std::vector<int> empty_row_indices;
  std::vector<int> empty_col_ptrs(1, 0);
  GatherAndBroadcastResult(rank, size, empty_values, empty_row_indices, empty_col_ptrs, total_cols, a.rows, result);
  return result;
}

bool KotelnikovaATaskALL::RunImpl() {
  const auto &[a, b] = GetInput();
  GetOutput() = MultiplyMatricesMPIOMP(a, b);
  return true;
}

bool KotelnikovaATaskALL::PostProcessingImpl() {
  return true;
}

}  // namespace kotelnikova_a_double_matr_mult
