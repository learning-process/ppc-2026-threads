#include "volkov_a_sparse_mat_mul_ccs/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cmath>
#include <tuple>
#include <utility>
#include <vector>

#include "volkov_a_sparse_mat_mul_ccs/common/include/common.hpp"

namespace volkov_a_sparse_mat_mul_ccs {

namespace {

template <typename Matrix>
void BroadcastMatrix(Matrix &mat, int root, MPI_Comm comm) {
  int rank = 0;
  MPI_Comm_rank(comm, &rank);

  MPI_Bcast(&mat.rows_count, 1, MPI_INT, root, comm);
  MPI_Bcast(&mat.cols_count, 1, MPI_INT, root, comm);
  MPI_Bcast(&mat.non_zeros, 1, MPI_INT, root, comm);

  if (rank != root) {
    mat.col_ptrs.resize(mat.cols_count + 1);
    mat.row_indices.resize(mat.non_zeros);
    mat.values.resize(mat.non_zeros);
  }

  MPI_Bcast(mat.col_ptrs.data(), mat.cols_count + 1, MPI_INT, root, comm);
  if (mat.non_zeros > 0) {
    MPI_Bcast(mat.row_indices.data(), mat.non_zeros, MPI_INT, root, comm);
    MPI_Bcast(mat.values.data(), mat.non_zeros, MPI_DOUBLE, root, comm);
  }
}

template <typename MatrixType>
void ProcessColumnOmp(int col_idx, const MatrixType &matrix_a, const MatrixType &matrix_b,
                      std::vector<double> &col_accumulator, std::vector<int> &local_row_indices,
                      std::vector<double> &local_values) {
  int b_start = matrix_b.col_ptrs[col_idx];
  int b_end = matrix_b.col_ptrs[col_idx + 1];

  for (int k = b_start; k < b_end; ++k) {
    int b_row = matrix_b.row_indices[k];
    double b_val = matrix_b.values[k];

    int a_start = matrix_a.col_ptrs[b_row];
    int a_end = matrix_a.col_ptrs[b_row + 1];

    for (int idx = a_start; idx < a_end; ++idx) {
      int a_row = matrix_a.row_indices[idx];
      double a_val = matrix_a.values[idx];
      col_accumulator[a_row] += a_val * b_val;
    }
  }

  for (int i = 0; i < matrix_a.rows_count; ++i) {
    if (std::abs(col_accumulator[i]) > 1e-10) {
      local_row_indices.push_back(i);
      local_values.push_back(col_accumulator[i]);
    }
    col_accumulator[i] = 0.0;
  }
}

void FlattenLocalData(int my_cols_count, const std::vector<std::vector<int>> &local_row_indices,
                      const std::vector<std::vector<double>> &local_values, std::vector<int> &my_col_sizes, int &my_nnz,
                      std::vector<int> &my_rows_flat, std::vector<double> &my_vals_flat) {
  my_nnz = 0;
  my_col_sizes.resize(my_cols_count);
  for (int i = 0; i < my_cols_count; ++i) {
    my_col_sizes[i] = static_cast<int>(local_row_indices[i].size());
    my_nnz += my_col_sizes[i];
  }

  my_rows_flat.resize(my_nnz);
  my_vals_flat.resize(my_nnz);
  int offset = 0;
  for (int i = 0; i < my_cols_count; ++i) {
    std::copy(local_row_indices[i].begin(), local_row_indices[i].end(), my_rows_flat.begin() + offset);
    std::copy(local_values[i].begin(), local_values[i].end(), my_vals_flat.begin() + offset);
    offset += my_col_sizes[i];
  }
}

void SendResultsToRoot(int my_cols_count, const std::vector<int> &my_col_sizes, int my_nnz,
                       const std::vector<int> &my_rows_flat, const std::vector<double> &my_vals_flat) {
  if (my_cols_count > 0) {
    MPI_Send(my_col_sizes.data(), my_cols_count, MPI_INT, 0, 0, MPI_COMM_WORLD);
    if (my_nnz > 0) {
      MPI_Send(my_rows_flat.data(), my_nnz, MPI_INT, 0, 1, MPI_COMM_WORLD);
      MPI_Send(my_vals_flat.data(), my_nnz, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
    }
  }
}

void ReceiveWorkerResults(int world_size, int chunk, int rem, std::vector<std::vector<int>> &all_row_indices,
                          std::vector<std::vector<double>> &all_values) {
  for (int proc = 1; proc < world_size; ++proc) {
    int proc_start_col = (proc * chunk) + std::min(proc, rem);
    int proc_end_col = proc_start_col + chunk + (proc < rem ? 1 : 0);
    int proc_cols_count = proc_end_col - proc_start_col;

    if (proc_cols_count <= 0) {
      continue;
    }

    std::vector<int> proc_col_sizes(proc_cols_count);
    MPI_Recv(proc_col_sizes.data(), proc_cols_count, MPI_INT, proc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    int proc_nnz = 0;
    for (int s : proc_col_sizes) {
      proc_nnz += s;
    }

    std::vector<int> proc_rows(proc_nnz);
    std::vector<double> proc_vals(proc_nnz);

    if (proc_nnz > 0) {
      MPI_Recv(proc_rows.data(), proc_nnz, MPI_INT, proc, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(proc_vals.data(), proc_nnz, MPI_DOUBLE, proc, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    int curr_offset = 0;
    for (int i = 0; i < proc_cols_count; ++i) {
      all_row_indices[proc_start_col + i].resize(proc_col_sizes[i]);
      all_values[proc_start_col + i].resize(proc_col_sizes[i]);

      std::copy(proc_rows.begin() + curr_offset, proc_rows.begin() + curr_offset + proc_col_sizes[i],
                all_row_indices[proc_start_col + i].begin());
      std::copy(proc_vals.begin() + curr_offset, proc_vals.begin() + curr_offset + proc_col_sizes[i],
                all_values[proc_start_col + i].begin());

      curr_offset += proc_col_sizes[i];
    }
  }
}

template <typename MatrixType>
void AssembleFinalMatrix(int total_cols, const std::vector<std::vector<int>> &all_row_indices,
                         const std::vector<std::vector<double>> &all_values, MatrixType &matrix_c) {
  matrix_c.col_ptrs.assign(total_cols + 1, 0);
  for (int j = 0; j < total_cols; ++j) {
    matrix_c.col_ptrs[j + 1] = matrix_c.col_ptrs[j] + static_cast<int>(all_row_indices[j].size());
  }

  matrix_c.non_zeros = matrix_c.col_ptrs.back();
  matrix_c.row_indices.resize(matrix_c.non_zeros);
  matrix_c.values.resize(matrix_c.non_zeros);

#pragma omp parallel for schedule(static) default(none) shared(total_cols, matrix_c, all_row_indices, all_values)
  for (int j = 0; j < total_cols; ++j) {
    int offset_c = matrix_c.col_ptrs[j];
    int size = static_cast<int>(all_row_indices[j].size());
    for (int k = 0; k < size; ++k) {
      matrix_c.row_indices[offset_c + k] = all_row_indices[j][k];
      matrix_c.values[offset_c + k] = all_values[j][k];
    }
  }
}

template <typename MatrixType>
void GatherResultsToRoot(int total_cols, int world_size, int chunk, int rem, int start_col, int my_cols_count,
                         std::vector<std::vector<int>> &local_row_indices,
                         std::vector<std::vector<double>> &local_values, MatrixType &matrix_c) {
  std::vector<std::vector<int>> all_row_indices(total_cols);
  std::vector<std::vector<double>> all_values(total_cols);

  for (int i = 0; i < my_cols_count; ++i) {
    all_row_indices[start_col + i] = std::move(local_row_indices[i]);
    all_values[start_col + i] = std::move(local_values[i]);
  }

  ReceiveWorkerResults(world_size, chunk, rem, all_row_indices, all_values);

  AssembleFinalMatrix(total_cols, all_row_indices, all_values, matrix_c);
}

}  // namespace

VolkovASparseMatMulCcsAll::VolkovASparseMatMulCcsAll(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool VolkovASparseMatMulCcsAll::ValidationImpl() {
  int world_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  if (world_rank == 0) {
    const auto &matrix_a = std::get<0>(GetInput());
    const auto &matrix_b = std::get<1>(GetInput());
    return (matrix_a.cols_count == matrix_b.rows_count);
  }
  return true;
}

bool VolkovASparseMatMulCcsAll::PreProcessingImpl() {
  return true;
}

bool VolkovASparseMatMulCcsAll::RunImpl() {
  int world_rank = 0;
  int world_size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  auto &matrix_a = std::get<0>(GetInput());
  auto &matrix_b = std::get<1>(GetInput());
  auto &matrix_c = GetOutput();

  BroadcastMatrix(matrix_a, 0, MPI_COMM_WORLD);
  BroadcastMatrix(matrix_b, 0, MPI_COMM_WORLD);

  if (world_rank == 0) {
    matrix_c.rows_count = matrix_a.rows_count;
    matrix_c.cols_count = matrix_b.cols_count;
  }

  int total_cols = matrix_b.cols_count;
  int chunk = total_cols / world_size;
  int rem = total_cols % world_size;
  int start_col = (world_rank * chunk) + std::min(world_rank, rem);
  int end_col = start_col + chunk + (world_rank < rem ? 1 : 0);
  int my_cols_count = end_col - start_col;

  std::vector<std::vector<int>> local_row_indices(my_cols_count);
  std::vector<std::vector<double>> local_values(my_cols_count);

#pragma omp parallel default(none) shared(matrix_a, matrix_b, start_col, end_col, local_row_indices, local_values)
  {
    std::vector<double> col_accumulator(matrix_a.rows_count, 0.0);
#pragma omp for schedule(dynamic)
    for (int j = start_col; j < end_col; ++j) {
      int local_j = j - start_col;
      ProcessColumnOmp(j, matrix_a, matrix_b, col_accumulator, local_row_indices[local_j], local_values[local_j]);
    }
  }

  if (world_rank == 0) {
    GatherResultsToRoot(total_cols, world_size, chunk, rem, start_col, my_cols_count, local_row_indices, local_values,
                        matrix_c);
  } else {
    int my_nnz = 0;
    std::vector<int> my_col_sizes;
    std::vector<int> my_rows_flat;
    std::vector<double> my_vals_flat;

    FlattenLocalData(my_cols_count, local_row_indices, local_values, my_col_sizes, my_nnz, my_rows_flat, my_vals_flat);
    SendResultsToRoot(my_cols_count, my_col_sizes, my_nnz, my_rows_flat, my_vals_flat);
  }

  return true;
}

bool VolkovASparseMatMulCcsAll::PostProcessingImpl() {
  auto &matrix_c = GetOutput();
  BroadcastMatrix(matrix_c, 0, MPI_COMM_WORLD);
  return true;
}

}  // namespace volkov_a_sparse_mat_mul_ccs
