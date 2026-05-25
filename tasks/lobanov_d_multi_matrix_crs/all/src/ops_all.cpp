#include "lobanov_d_multi_matrix_crs/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

#include "lobanov_d_multi_matrix_crs/common/include/common.hpp"
#include "util/include/util.hpp"

namespace lobanov_d_multi_matrix_crs {

void LobanovMultyMatrixALL::SortIndices(std::vector<int> &vec) {
  std::ranges::sort(vec);
}

LobanovMultyMatrixALL::LobanovMultyMatrixALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    GetInput() = in;
  }
}

bool LobanovMultyMatrixALL::ValidationImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int valid_flag = 0;
  if (rank == 0) {
    const auto &input = GetInput();
    const auto &mat_a = input.first;
    const auto &mat_b = input.second;
    if (mat_a.column_count == mat_b.row_count && mat_a.row_count > 0 && mat_b.column_count > 0) {
      valid_flag = 1;
    }
  }
  MPI_Bcast(&valid_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
  return valid_flag == 1;
}

bool LobanovMultyMatrixALL::PreProcessingImpl() {
  return true;
}

void LobanovMultyMatrixALL::DistributeSparseMatrix(CompressedRowMatrix &mat, int root, int rows, int cols) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::array<int, 2> sizes = {0, 0};
  if (rank == root) {
    sizes[0] = static_cast<int>(mat.row_pointer_data.size());
    sizes[1] = static_cast<int>(mat.value_data.size());
  }
  MPI_Bcast(sizes.data(), 2, MPI_INT, root, MPI_COMM_WORLD);
  if (rank != root) {
    mat.row_pointer_data.resize(sizes[0]);
    mat.value_data.resize(sizes[1]);
    mat.column_index_data.resize(sizes[1]);
  }
  if (sizes[0] > 0) {
    MPI_Bcast(mat.row_pointer_data.data(), sizes[0], MPI_INT, root, MPI_COMM_WORLD);
  }
  if (sizes[1] > 0) {
    MPI_Bcast(mat.value_data.data(), sizes[1], MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Bcast(mat.column_index_data.data(), sizes[1], MPI_INT, root, MPI_COMM_WORLD);
  }
  mat.row_count = rows;
  mat.column_count = cols;
  mat.non_zero_count = static_cast<int>(mat.value_data.size());
}

CompressedRowMatrix LobanovMultyMatrixALL::TransposeSparseMatrix(const CompressedRowMatrix &src) {
  CompressedRowMatrix dst;
  dst.row_count = src.column_count;
  dst.column_count = src.row_count;
  dst.row_pointer_data.assign(dst.row_count + 1, 0);

  for (int col : src.column_index_data) {
    dst.row_pointer_data[col + 1]++;
  }
  for (int i = 0; i < dst.row_count; ++i) {
    dst.row_pointer_data[i + 1] += dst.row_pointer_data[i];
  }
  dst.value_data.resize(src.value_data.size());
  dst.column_index_data.resize(src.column_index_data.size());
  dst.non_zero_count = src.non_zero_count;

  std::vector<int> cursor = dst.row_pointer_data;
  for (int i = 0; i < src.row_count; ++i) {
    for (int j = src.row_pointer_data[i]; j < src.row_pointer_data[i + 1]; ++j) {
      int col = src.column_index_data[j];
      int pos = cursor[col];
      cursor[col]++;
      dst.value_data[pos] = src.value_data[j];
      dst.column_index_data[pos] = i;
    }
  }
  return dst;
}

void LobanovMultyMatrixALL::ComputeLocalProduct(const CompressedRowMatrix &a, const CompressedRowMatrix &b_tr,
                                                int start_row, int local_rows, std::vector<int> &row_nnz_counts,
                                                std::vector<double> &packed_vals, std::vector<int> &packed_cols) {
  if (local_rows <= 0) {
    return;
  }
  int result_cols = b_tr.row_count;

  std::vector<std::vector<double>> row_vals(local_rows);
  std::vector<std::vector<int>> row_cols(local_rows);

#pragma omp parallel default(none) shared(a, b_tr, start_row, local_rows, row_vals, row_cols, row_nnz_counts, \
                                              result_cols) num_threads(ppc::util::GetNumThreads())
  {
    std::vector<int> marker(result_cols, -1);
    std::vector<double> accumulator(result_cols, 0.0);
    std::vector<int> active_columns;

#pragma omp for schedule(dynamic)
    for (int i = 0; i < local_rows; ++i) {
      int global_row = start_row + i;
      active_columns.clear();

      for (int idx = a.row_pointer_data[global_row]; idx < a.row_pointer_data[global_row + 1]; ++idx) {
        int col_a = a.column_index_data[idx];
        double val_a = a.value_data[idx];

        for (int j = b_tr.row_pointer_data[col_a]; j < b_tr.row_pointer_data[col_a + 1]; ++j) {
          int col_res = b_tr.column_index_data[j];
          double contrib = val_a * b_tr.value_data[j];
          if (marker[col_res] != i) {
            marker[col_res] = i;
            active_columns.push_back(col_res);
            accumulator[col_res] = contrib;
          } else {
            accumulator[col_res] += contrib;
          }
        }
      }

      row_nnz_counts[i] = static_cast<int>(active_columns.size());
      SortIndices(active_columns);
      for (int col : active_columns) {
        row_vals[i].push_back(accumulator[col]);
        row_cols[i].push_back(col);
        accumulator[col] = 0.0;
      }
    }
  }

  for (int i = 0; i < local_rows; ++i) {
    packed_vals.insert(packed_vals.end(), row_vals[i].begin(), row_vals[i].end());
    packed_cols.insert(packed_cols.end(), row_cols[i].begin(), row_cols[i].end());
  }
}

void LobanovMultyMatrixALL::MergeLocalResults(int rank, int comm_size, int total_rows, int result_cols, int local_rows,
                                              CompressedRowMatrix &result_mat, const std::vector<int> &row_nnz_counts,
                                              const std::vector<double> &packed_vals,
                                              const std::vector<int> &packed_cols) {
  int local_total_nnz = static_cast<int>(packed_vals.size());
  std::vector<int> global_nnz_counts(comm_size);
  std::vector<int> global_row_counts(comm_size);

  MPI_Gather(&local_total_nnz, 1, MPI_INT, global_nnz_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&local_rows, 1, MPI_INT, global_row_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

  // Создаем неконстантные копии для MPI (mpi не меняет данные)
  std::vector<double> packed_vals_copy = packed_vals;
  std::vector<int> packed_cols_copy = packed_cols;
  std::vector<int> row_nnz_counts_copy = row_nnz_counts;

  if (rank == 0) {
    result_mat.row_count = total_rows;
    result_mat.column_count = result_cols;

    std::vector<int> nnz_offsets(comm_size, 0);
    std::vector<int> row_offsets(comm_size, 0);
    int total_nnz = 0;
    for (int proc = 0; proc < comm_size; ++proc) {
      nnz_offsets[proc] = total_nnz;
      total_nnz += global_nnz_counts[proc];
      if (proc > 0) {
        row_offsets[proc] = row_offsets[proc - 1] + global_row_counts[proc - 1];
      }
    }

    result_mat.value_data.resize(total_nnz);
    result_mat.column_index_data.resize(total_nnz);
    result_mat.row_pointer_data.assign(static_cast<size_t>(total_rows) + 1, 0);
    result_mat.non_zero_count = total_nnz;

    MPI_Gatherv(packed_vals_copy.data(), local_total_nnz, MPI_DOUBLE, result_mat.value_data.data(),
                global_nnz_counts.data(), nnz_offsets.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gatherv(packed_cols_copy.data(), local_total_nnz, MPI_INT, result_mat.column_index_data.data(),
                global_nnz_counts.data(), nnz_offsets.data(), MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<int> all_row_nnz(total_rows);
    MPI_Gatherv(row_nnz_counts_copy.data(), local_rows, MPI_INT, all_row_nnz.data(), global_row_counts.data(),
                row_offsets.data(), MPI_INT, 0, MPI_COMM_WORLD);
    for (int i = 0; i < total_rows; ++i) {
      result_mat.row_pointer_data[i + 1] = result_mat.row_pointer_data[i] + all_row_nnz[i];
    }
  } else {
    MPI_Gatherv(packed_vals_copy.data(), local_total_nnz, MPI_DOUBLE, nullptr, nullptr, nullptr, MPI_DOUBLE, 0,
                MPI_COMM_WORLD);
    MPI_Gatherv(packed_cols_copy.data(), local_total_nnz, MPI_INT, nullptr, nullptr, nullptr, MPI_INT, 0,
                MPI_COMM_WORLD);
    MPI_Gatherv(row_nnz_counts_copy.data(), local_rows, MPI_INT, nullptr, nullptr, nullptr, MPI_INT, 0, MPI_COMM_WORLD);
  }

  DistributeSparseMatrix(result_mat, 0, total_rows, result_cols);
}

bool LobanovMultyMatrixALL::RunImpl() {
  int comm_size = 0;
  int rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::array<int, 3> dimensions = {0, 0, 0};
  if (rank == 0) {
    const auto &mat_a = GetInput().first;
    const auto &mat_b = GetInput().second;
    dimensions[0] = mat_a.row_count;
    dimensions[1] = mat_a.column_count;
    dimensions[2] = mat_b.column_count;
  }
  MPI_Bcast(dimensions.data(), 3, MPI_INT, 0, MPI_COMM_WORLD);
  int a_rows = dimensions[0];
  int a_cols = dimensions[1];
  int b_cols = dimensions[2];

  CompressedRowMatrix matrix_a;
  CompressedRowMatrix matrix_b;
  if (rank == 0) {
    matrix_a = GetInput().first;
    matrix_b = GetInput().second;
  }
  DistributeSparseMatrix(matrix_a, 0, a_rows, a_cols);
  DistributeSparseMatrix(matrix_b, 0, a_cols, b_cols);

  CompressedRowMatrix matrix_b_transposed = TransposeSparseMatrix(matrix_b);

  int base_chunk = a_rows / comm_size;
  int remainder = a_rows % comm_size;
  int start_row = (rank * base_chunk) + std::min(rank, remainder);
  int local_rows = base_chunk + (rank < remainder ? 1 : 0);

  std::vector<int> local_row_nnz(local_rows, 0);
  std::vector<double> flat_values;
  std::vector<int> flat_columns;

  ComputeLocalProduct(matrix_a, matrix_b_transposed, start_row, local_rows, local_row_nnz, flat_values, flat_columns);

  CompressedRowMatrix result_matrix;
  MergeLocalResults(rank, comm_size, a_rows, b_cols, local_rows, result_matrix, local_row_nnz, flat_values,
                    flat_columns);

  if (rank == 0) {
    GetOutput() = result_matrix;
  }
  return true;
}

bool LobanovMultyMatrixALL::PostProcessingImpl() {
  return true;
}

}  // namespace lobanov_d_multi_matrix_crs
