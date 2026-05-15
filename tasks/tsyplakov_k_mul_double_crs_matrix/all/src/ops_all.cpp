#include "tsyplakov_k_mul_double_crs_matrix/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

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

  for (int idx_a = a.row_ptr[row]; idx_a < a.row_ptr[row + 1]; ++idx_a) {
    const int k = a.col_index[idx_a];
    const double val_a = a.values[idx_a];

    for (int idx_b = b.row_ptr[k]; idx_b < b.row_ptr[k + 1]; ++idx_b) {
      const int j = b.col_index[idx_b];
      acc[j] += val_a * b.values[idx_b];
    }
  }

  values.reserve(acc.size());
  cols.reserve(acc.size());

  for (const auto &[col, val] : acc) {
    if (std::fabs(val) > 1e-12) {
      cols.push_back(col);
      values.push_back(val);
    }
  }
}

}  // namespace

TsyplakovKTestTaskALL::TsyplakovKTestTaskALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool TsyplakovKTestTaskALL::ValidationImpl() {
  const auto &input = GetInput();
  return input.a.cols == input.b.rows;
}

bool TsyplakovKTestTaskALL::PreProcessingImpl() {
  return true;
}

bool TsyplakovKTestTaskALL::RunImpl() {
  const auto &input = GetInput();

  const auto &a = input.a;
  const auto &b = input.b;

  int rank = 0;
  int size = 1;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const int rows = a.rows;

  const int rows_per_proc = rows / size;
  const int remainder = rows % size;

  const int start_row = (rank * rows_per_proc) + std::min(rank, remainder);

  const int local_rows = rows_per_proc + (rank < remainder ? 1 : 0);

  const int end_row = start_row + local_rows;

  std::vector<std::vector<double>> local_values(local_rows);
  std::vector<std::vector<int>> local_cols(local_rows);

#pragma omp parallel for default(none) shared(a, b, local_values, local_cols, start_row, end_row)
  for (int row = start_row; row < end_row; ++row) {
    const int local_idx = row - start_row;

    ComputeRow(a, b, row, local_values[local_idx], local_cols[local_idx]);
  }

  SparseMatrixCRS local_matrix(local_rows, b.cols);

  for (int i = 0; i < local_rows; ++i) {
    local_matrix.row_ptr[i + 1] = local_matrix.row_ptr[i] + static_cast<int>(local_values[i].size());
  }

  const int local_nnz = local_matrix.row_ptr[local_rows];

  local_matrix.values.reserve(local_nnz);
  local_matrix.col_index.reserve(local_nnz);

  for (int i = 0; i < local_rows; ++i) {
    local_matrix.values.insert(local_matrix.values.end(), local_values[i].begin(), local_values[i].end());

    local_matrix.col_index.insert(local_matrix.col_index.end(), local_cols[i].begin(), local_cols[i].end());
  }

  std::vector<int> recv_nnz(size);
  std::vector<int> recv_rows(size);

  MPI_Gather(&local_nnz, 1, MPI_INT, recv_nnz.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Gather(&local_rows, 1, MPI_INT, recv_rows.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<int> displs_nnz(size, 0);
  std::vector<int> displs_rows(size, 0);

  int total_nnz = 0;

  if (rank == 0) {
    for (int i = 1; i < size; ++i) {
      displs_nnz[i] = displs_nnz[i - 1] + recv_nnz[i - 1];

      displs_rows[i] = displs_rows[i - 1] + recv_rows[i - 1];
    }

    total_nnz = displs_nnz[size - 1] + recv_nnz[size - 1];
  }

  SparseMatrixCRS result_matrix;

  double *recv_values = nullptr;
  int *recv_cols = nullptr;
  int *gathered_rows_ptr = nullptr;

  std::vector<int> gathered_row_sizes;

  if (rank == 0) {
    result_matrix.rows = rows;
    result_matrix.cols = b.cols;

    result_matrix.values.resize(total_nnz);
    result_matrix.col_index.resize(total_nnz);
    result_matrix.row_ptr.resize(rows + 1, 0);

    recv_values = result_matrix.values.data();
    recv_cols = result_matrix.col_index.data();

    gathered_row_sizes.resize(rows);
    gathered_rows_ptr = gathered_row_sizes.data();
  }

  MPI_Gatherv(local_matrix.values.data(), local_nnz, MPI_DOUBLE, recv_values, recv_nnz.data(), displs_nnz.data(),
              MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Gatherv(local_matrix.col_index.data(), local_nnz, MPI_INT, recv_cols, recv_nnz.data(), displs_nnz.data(), MPI_INT,
              0, MPI_COMM_WORLD);

  std::vector<int> local_row_sizes(local_rows);

  for (int i = 0; i < local_rows; ++i) {
    local_row_sizes[i] = static_cast<int>(local_values[i].size());
  }

  MPI_Gatherv(local_row_sizes.data(), local_rows, MPI_INT, gathered_rows_ptr, recv_rows.data(), displs_rows.data(),
              MPI_INT, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    for (int i = 0; i < rows; ++i) {
      result_matrix.row_ptr[i + 1] = result_matrix.row_ptr[i] + gathered_row_sizes[i];
    }

    GetOutput() = std::move(result_matrix);
  }

  return true;
}

bool TsyplakovKTestTaskALL::PostProcessingImpl() {
  return true;
}

}  // namespace tsyplakov_k_mul_double_crs_matrix
