#include "tsyplakov_k_mul_double_crs_matrix/all/include/ops_all.hpp"

#include <mpi.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <utility>
#include <vector>

#include "tsyplakov_k_mul_double_crs_matrix/common/include/common.hpp"

namespace tsyplakov_k_mul_double_crs_matrix {

namespace {

void ComputeRow(const SparseMatrixCRS &a, const SparseMatrixCRS &b, int row, std::vector<double> &values,
                std::vector<int> &cols) {
  std::unordered_map<int, double> acc;

  for (int ia = a.row_ptr[row]; ia < a.row_ptr[row + 1]; ++ia) {
    const int k = a.col_index[ia];
    const double va = a.values[ia];

    for (int ib = b.row_ptr[k]; ib < b.row_ptr[k + 1]; ++ib) {
      acc[b.col_index[ib]] += va * b.values[ib];
    }
  }

  for (const auto &[c, v] : acc) {
    if (std::fabs(v) > 1e-12) {
      cols.push_back(c);
      values.push_back(v);
    }
  }
}

}  // namespace

TsyplakovKTestTaskALL::TsyplakovKTestTaskALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool TsyplakovKTestTaskALL::ValidationImpl() {
  const auto &in = GetInput();
  return in.a.cols == in.b.rows;
}

bool TsyplakovKTestTaskALL::PreProcessingImpl() {
  return true;
}

bool TsyplakovKTestTaskALL::RunImpl() {
  const auto &a = GetInput().a;
  const auto &b = GetInput().b;

  int rank = 0;
  int size = 1;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const int n = a.rows;

  const int base = n / size;
  const int rem = n % size;

  const int start = (rank * base) + std::min(rank, rem);
  const int local_n = base + (rank < rem ? 1 : 0);

  std::vector<std::vector<double>> loc_vals(local_n);
  std::vector<std::vector<int>> loc_cols(local_n);

  tbb::parallel_for(tbb::blocked_range<int>(0, local_n), [&](const tbb::blocked_range<int> &r) {
    for (int i = r.begin(); i < r.end(); ++i) {
      ComputeRow(a, b, start + i, loc_vals[i], loc_cols[i]);
    }
  });

  std::vector<int> local_sizes(local_n);
  for (int i = 0; i < local_n; ++i) {
    local_sizes[i] = static_cast<int>(loc_vals[i].size());
  }

  std::vector<int> rows_per_proc(size);
  MPI_Gather(&local_n, 1, MPI_INT, rows_per_proc.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<int> row_displs(size, 0);
  if (rank == 0) {
    for (int i = 1; i < size; ++i) {
      row_displs[i] = row_displs[i - 1] + rows_per_proc[i - 1];
    }
  }

  std::vector<int> flat_sizes(n, 0);

  MPI_Gatherv(local_sizes.data(), local_n, MPI_INT, flat_sizes.data(), rows_per_proc.data(), row_displs.data(), MPI_INT,
              0, MPI_COMM_WORLD);

  MPI_Bcast(flat_sizes.data(), n, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<double> local_v;
  std::vector<int> local_c;

  int local_nnz = 0;
  for (int i = 0; i < local_n; ++i) {
    local_nnz += static_cast<int>(loc_vals[i].size());
    local_v.insert(local_v.end(), loc_vals[i].begin(), loc_vals[i].end());
    local_c.insert(local_c.end(), loc_cols[i].begin(), loc_cols[i].end());
  }

  std::vector<int> nnz_counts(size);
  MPI_Gather(&local_nnz, 1, MPI_INT, nnz_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<int> nnz_disp(size, 0);
  int total_nnz = 0;
  if (rank == 0) {
    for (int i = 1; i < size; ++i) {
      nnz_disp[i] = nnz_disp[i - 1] + nnz_counts[i - 1];
    }
    total_nnz = nnz_disp[size - 1] + nnz_counts[size - 1];
  }

  MPI_Bcast(nnz_counts.data(), size, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(nnz_disp.data(), size, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&total_nnz, 1, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<double> out_v(total_nnz);
  std::vector<int> out_c(total_nnz);
  std::vector<int> out_row_ptr(n + 1, 0);

  MPI_Allgatherv(local_v.data(), local_nnz, MPI_DOUBLE, out_v.data(), nnz_counts.data(), nnz_disp.data(), MPI_DOUBLE,
                 MPI_COMM_WORLD);

  MPI_Allgatherv(local_c.data(), local_nnz, MPI_INT, out_c.data(), nnz_counts.data(), nnz_disp.data(), MPI_INT,
                 MPI_COMM_WORLD);

  for (int i = 0; i < n; ++i) {
    out_row_ptr[i + 1] = out_row_ptr[i] + flat_sizes[i];
  }

  SparseMatrixCRS res(n, b.cols);
  res.values = std::move(out_v);
  res.col_index = std::move(out_c);
  res.row_ptr = std::move(out_row_ptr);
  GetOutput() = std::move(res);

  return true;
}

bool TsyplakovKTestTaskALL::PostProcessingImpl() {
  return true;
}

}  // namespace tsyplakov_k_mul_double_crs_matrix
