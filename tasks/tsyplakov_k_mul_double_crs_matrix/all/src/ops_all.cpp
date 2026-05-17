#include "tsyplakov_k_mul_double_crs_matrix/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include "tsyplakov_k_mul_double_crs_matrix/common/include/common.hpp"

namespace tsyplakov_k_mul_double_crs_matrix {

namespace {

void ComputeRow(const SparseMatrixCRS& a,
                const SparseMatrixCRS& b,
                int row,
                std::vector<double>& values,
                std::vector<int>& cols) {
  std::vector<double> accumulator(b.cols, 0.0);
  std::vector<int> used_cols;

  for (int idx_a = a.row_ptr[row]; idx_a < a.row_ptr[row + 1]; ++idx_a) {
    const int k = a.col_index[idx_a];
    const double val_a = a.values[idx_a];

    for (int idx_b = b.row_ptr[k]; idx_b < b.row_ptr[k + 1]; ++idx_b) {
      const int col = b.col_index[idx_b];

      if (std::fabs(accumulator[col]) < 1e-12) {
        used_cols.push_back(col);
      }

      accumulator[col] += val_a * b.values[idx_b];
    }
  }

  values.reserve(used_cols.size());
  cols.reserve(used_cols.size());

  for (const int col : used_cols) {
    if (std::fabs(accumulator[col]) > 1e-12) {
      cols.push_back(col);
      values.push_back(accumulator[col]);
    }
  }
}

}  // namespace

TsyplakovKTestTaskALL::TsyplakovKTestTaskALL(const InType& in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool TsyplakovKTestTaskALL::ValidationImpl() {
  const auto& input = GetInput();

  return input.a.cols == input.b.rows;
}

bool TsyplakovKTestTaskALL::PreProcessingImpl() {
  return true;
}

bool TsyplakovKTestTaskALL::RunImpl() {
  const auto& input = GetInput();

  const auto& a = input.a;
  const auto& b = input.b;

  int rank = 0;
  int size = 1;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const int rows = a.rows;

  const int rows_per_process = rows / size;
  const int remainder = rows % size;

  const int start_row =
      (rank * rows_per_process) + std::min(rank, remainder);

  const int local_rows =
      rows_per_process + ((rank < remainder) ? 1 : 0);

  const int end_row = start_row + local_rows;

  std::vector<std::vector<double>> local_values(local_rows);
  std::vector<std::vector<int>> local_cols(local_rows);

  tbb::parallel_for(
      tbb::blocked_range<int>(start_row, end_row),
      [&](const tbb::blocked_range<int>& range) {
        for (int row = range.begin(); row < range.end(); ++row) {
          const int local_index = row - start_row;

          ComputeRow(a,
                     b,
                     row,
                     local_values[local_index],
                     local_cols[local_index]);
        }
      });

  SparseMatrixCRS local_matrix(local_rows, b.cols);

  for (int i = 0; i < local_rows; ++i) {
    local_matrix.row_ptr[i + 1] =
        local_matrix.row_ptr[i] +
        static_cast<int>(local_values[i].size());
  }

  const int local_nnz = local_matrix.row_ptr[local_rows];

  local_matrix.values.reserve(local_nnz);
  local_matrix.col_index.reserve(local_nnz);

  for (int i = 0; i < local_rows; ++i) {
    local_matrix.values.insert(local_matrix.values.end(),
                               local_values[i].begin(),
                               local_values[i].end());

    local_matrix.col_index.insert(local_matrix.col_index.end(),
                                  local_cols[i].begin(),
                                  local_cols[i].end());
  }

  std::vector<int> recv_nnz(size);
  std::vector<int> recv_rows(size);

  MPI_Gather(&local_nnz,
             1,
             MPI_INT,
             recv_nnz.data(),
             1,
             MPI_INT,
             0,
             MPI_COMM_WORLD);

  MPI_Gather(&local_rows,
             1,
             MPI_INT,
             recv_rows.data(),
             1,
             MPI_INT,
             0,
             MPI_COMM_WORLD);

  std::vector<int> displs_nnz(size, 0);
  std::vector<int> displs_rows(size, 0);

  int total_nnz = 0;

  if (rank == 0) {
    for (int i = 1; i < size; ++i) {
      displs_nnz[i] =
          displs_nnz[i - 1] + recv_nnz[i - 1];

      displs_rows[i] =
          displs_rows[i - 1] + recv_rows[i - 1];
    }

    total_nnz =
        displs_nnz[size - 1] + recv_nnz[size - 1];
  }

  SparseMatrixCRS result_matrix;

  double* recv_values = nullptr;
  int* recv_cols = nullptr;
  int* recv_row_sizes = nullptr;

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
    recv_row_sizes = gathered_row_sizes.data();
  }

  double* send_values =
      local_matrix.values.empty()
          ? nullptr
          : local_matrix.values.data();

  int* send_cols =
      local_matrix.col_index.empty()
          ? nullptr
          : local_matrix.col_index.data();

  MPI_Gatherv(send_values,
              local_nnz,
              MPI_DOUBLE,
              recv_values,
              recv_nnz.data(),
              displs_nnz.data(),
              MPI_DOUBLE,
              0,
              MPI_COMM_WORLD);

  MPI_Gatherv(send_cols,
              local_nnz,
              MPI_INT,
              recv_cols,
              recv_nnz.data(),
              displs_nnz.data(),
              MPI_INT,
              0,
              MPI_COMM_WORLD);

  std::vector<int> local_row_sizes(local_rows);

  for (int i = 0; i < local_rows; ++i) {
    local_row_sizes[i] =
        static_cast<int>(local_values[i].size());
  }

  int* send_row_sizes =
      local_row_sizes.empty()
          ? nullptr
          : local_row_sizes.data();

  MPI_Gatherv(send_row_sizes,
              local_rows,
              MPI_INT,
              recv_row_sizes,
              recv_rows.data(),
              displs_rows.data(),
              MPI_INT,
              0,
              MPI_COMM_WORLD);

  if (rank == 0) {
    for (int i = 0; i < rows; ++i) {
      result_matrix.row_ptr[i + 1] =
          result_matrix.row_ptr[i] +
          gathered_row_sizes[i];
    }
  }

  int global_nnz = 0;

  if (rank == 0) {
    global_nnz =
        static_cast<int>(result_matrix.values.size());
  }

  MPI_Bcast(&global_nnz,
            1,
            MPI_INT,
            0,
            MPI_COMM_WORLD);

  if (rank != 0) {
    result_matrix.rows = rows;
    result_matrix.cols = b.cols;

    result_matrix.values.resize(global_nnz);
    result_matrix.col_index.resize(global_nnz);
    result_matrix.row_ptr.resize(rows + 1);
  }

  MPI_Bcast(result_matrix.row_ptr.data(),
            rows + 1,
            MPI_INT,
            0,
            MPI_COMM_WORLD);

  if (global_nnz > 0) {
    MPI_Bcast(result_matrix.values.data(),
              global_nnz,
              MPI_DOUBLE,
              0,
              MPI_COMM_WORLD);

    MPI_Bcast(result_matrix.col_index.data(),
              global_nnz,
              MPI_INT,
              0,
              MPI_COMM_WORLD);
  }

  GetOutput() = std::move(result_matrix);

  return true;
}

bool TsyplakovKTestTaskALL::PostProcessingImpl() {
  return true;
}

}  // namespace tsyplakov_k_mul_double_crs_matrix