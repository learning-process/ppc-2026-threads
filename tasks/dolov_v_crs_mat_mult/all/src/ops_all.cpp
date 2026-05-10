#include "dolov_v_crs_mat_mult/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

#include "dolov_v_crs_mat_mult/common/include/common.hpp"
#include "util/include/util.hpp"

namespace dolov_v_crs_mat_mult {

DolovVCrsMatMultAll::DolovVCrsMatMultAll(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    GetInput() = in;
  }
}

bool DolovVCrsMatMultAll::ValidationImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    const auto &input_data = GetInput();
    if (input_data.size() != 2) {
      return false;
    }
    const auto &matrix_a = input_data[0];
    const auto &matrix_b = input_data[1];
    return matrix_a.num_cols == matrix_b.num_rows && matrix_a.num_rows > 0 && matrix_b.num_cols > 0;
  }
  return true;
}

bool DolovVCrsMatMultAll::PreProcessingImpl() {
  return true;
}

SparseMatrix DolovVCrsMatMultAll::TransposeMatrix(const SparseMatrix &matrix) {
  SparseMatrix transposed;
  transposed.num_rows = matrix.num_cols;
  transposed.num_cols = matrix.num_rows;
  transposed.row_pointers.assign(transposed.num_rows + 1, 0);

  for (int col_idx : matrix.col_indices) {
    transposed.row_pointers[col_idx + 1]++;
  }
  for (int i = 0; i < transposed.num_rows; ++i) {
    transposed.row_pointers[i + 1] += transposed.row_pointers[i];
  }

  transposed.values.resize(matrix.values.size());
  transposed.col_indices.resize(matrix.col_indices.size());

  std::vector<int> current_pos = transposed.row_pointers;
  for (int i = 0; i < matrix.num_rows; ++i) {
    for (int j = matrix.row_pointers[i]; j < matrix.row_pointers[i + 1]; ++j) {
      int col = matrix.col_indices[j];
      int dest_idx = current_pos[col]++;
      transposed.values[dest_idx] = matrix.values[j];
      transposed.col_indices[dest_idx] = i;
    }
  }
  return transposed;
}

double DolovVCrsMatMultAll::DotProduct(const SparseMatrix &matrix_a, int row_a, const SparseMatrix &matrix_b_t,
                                       int row_b) {
  double sum = 0.0;
  int ptr_a = matrix_a.row_pointers[row_a];
  int ptr_b = matrix_b_t.row_pointers[row_b];
  const int end_a = matrix_a.row_pointers[row_a + 1];
  const int end_b = matrix_b_t.row_pointers[row_b + 1];

  while (ptr_a < end_a && ptr_b < end_b) {
    if (matrix_a.col_indices[ptr_a] == matrix_b_t.col_indices[ptr_b]) {
      sum += matrix_a.values[ptr_a] * matrix_b_t.values[ptr_b];
      ptr_a++;
      ptr_b++;
    } else if (matrix_a.col_indices[ptr_a] < matrix_b_t.col_indices[ptr_b]) {
      ptr_a++;
    } else {
      ptr_b++;
    }
  }
  return sum;
}

bool DolovVCrsMatMultAll::RunImpl() {
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  SparseMatrix local_a;
  SparseMatrix local_b_t;
  int sizes[4] = {0, 0, 0, 0};
  int a_nnz = 0;
  int b_t_nnz = 0;

  if (rank == 0) {
    const auto &matrix_a = GetInput()[0];
    const auto &matrix_b = GetInput()[1];

    local_a = matrix_a;
    local_b_t = TransposeMatrix(matrix_b);

    sizes[0] = local_a.num_rows;
    sizes[1] = local_a.num_cols;
    sizes[2] = local_b_t.num_rows;
    sizes[3] = local_b_t.num_cols;
    a_nnz = static_cast<int>(local_a.values.size());
    b_t_nnz = static_cast<int>(local_b_t.values.size());
  }

  MPI_Bcast(sizes, 4, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&a_nnz, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&b_t_nnz, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    local_a.num_rows = sizes[0];
    local_a.num_cols = sizes[1];
    local_a.row_pointers.resize(sizes[0] + 1);
    local_a.col_indices.resize(a_nnz);
    local_a.values.resize(a_nnz);

    local_b_t.num_rows = sizes[2];
    local_b_t.num_cols = sizes[3];
    local_b_t.row_pointers.resize(sizes[2] + 1);
    local_b_t.col_indices.resize(b_t_nnz);
    local_b_t.values.resize(b_t_nnz);
  }

  MPI_Bcast(local_a.row_pointers.data(), sizes[0] + 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (a_nnz > 0) {
    MPI_Bcast(local_a.col_indices.data(), a_nnz, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(local_a.values.data(), a_nnz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }

  MPI_Bcast(local_b_t.row_pointers.data(), sizes[2] + 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (b_t_nnz > 0) {
    MPI_Bcast(local_b_t.col_indices.data(), b_t_nnz, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(local_b_t.values.data(), b_t_nnz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }

  int a_rows = sizes[0];
  int chunk = a_rows / size;
  int remainder = a_rows % size;
  int local_start = rank * chunk + std::min(rank, remainder);
  int local_rows = chunk + (rank < remainder ? 1 : 0);

  std::vector<std::vector<double>> temp_values(local_rows);
  std::vector<std::vector<int>> temp_cols(local_rows);
  std::vector<int> local_nnz_per_row(local_rows, 0);

  int num_threads = ppc::util::GetNumThreads();
  if (num_threads <= 0) {
    num_threads = 1;
  }

#pragma omp parallel for default(none)                                                                    \
    shared(local_a, local_b_t, temp_values, temp_cols, local_nnz_per_row, local_rows, local_start, sizes) \
    num_threads(num_threads) schedule(dynamic)
  for (int i = 0; i < local_rows; ++i) {
    int global_row = local_start + i;
    std::vector<double> thread_vals;
    std::vector<int> thread_cols;

    for (int j = 0; j < sizes[2]; ++j) {
      double sum = DolovVCrsMatMultAll::DotProduct(local_a, global_row, local_b_t, j);
      if (std::abs(sum) > 1e-15) {
        thread_vals.push_back(sum);
        thread_cols.push_back(j);
      }
    }
    temp_values[i] = std::move(thread_vals);
    temp_cols[i] = std::move(thread_cols);
    local_nnz_per_row[i] = static_cast<int>(temp_values[i].size());
  }

  std::vector<double> flat_vals;
  std::vector<int> flat_cols;
  int local_nnz_total = 0;

  for (int i = 0; i < local_rows; ++i) {
    local_nnz_total += local_nnz_per_row[i];
  }

  flat_vals.reserve(local_nnz_total);
  flat_cols.reserve(local_nnz_total);
  for (int i = 0; i < local_rows; ++i) {
    flat_vals.insert(flat_vals.end(), temp_values[i].begin(), temp_values[i].end());
    flat_cols.insert(flat_cols.end(), temp_cols[i].begin(), temp_cols[i].end());
  }

  std::vector<int> recv_nnz_counts;
  std::vector<int> recv_nnz_displs;
  std::vector<int> rows_counts;
  std::vector<int> rows_displs;

  if (rank == 0) {
    recv_nnz_counts.resize(size);
    recv_nnz_displs.resize(size, 0);
    rows_counts.resize(size);
    rows_displs.resize(size, 0);

    for (int p = 0; p < size; ++p) {
      rows_counts[p] = chunk + (p < remainder ? 1 : 0);
      rows_displs[p] = p * chunk + std::min(p, remainder);
    }
  }

  MPI_Gather(&local_nnz_total, 1, MPI_INT, recv_nnz_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

  int total_global_nnz = 0;
  if (rank == 0) {
    for (int p = 0; p < size; ++p) {
      recv_nnz_displs[p] = total_global_nnz;
      total_global_nnz += recv_nnz_counts[p];
    }
  }

  std::vector<int> global_nnz_per_row;
  if (rank == 0) {
    global_nnz_per_row.resize(a_rows);
  }

  MPI_Gatherv(local_nnz_per_row.data(), local_rows, MPI_INT, global_nnz_per_row.data(), rows_counts.data(),
              rows_displs.data(), MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<double> global_vals;
  std::vector<int> global_cols;
  if (rank == 0) {
    global_vals.resize(total_global_nnz);
    global_cols.resize(total_global_nnz);
  }

  MPI_Gatherv(flat_vals.data(), local_nnz_total, MPI_DOUBLE, global_vals.data(), recv_nnz_counts.data(),
              recv_nnz_displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Gatherv(flat_cols.data(), local_nnz_total, MPI_INT, global_cols.data(), recv_nnz_counts.data(),
              recv_nnz_displs.data(), MPI_INT, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    SparseMatrix res;
    res.num_rows = a_rows;
    res.num_cols = sizes[2];
    res.row_pointers.assign(a_rows + 1, 0);

    for (int i = 0; i < a_rows; ++i) {
      res.row_pointers[i + 1] = res.row_pointers[i] + global_nnz_per_row[i];
    }
    res.values = std::move(global_vals);
    res.col_indices = std::move(global_cols);

    GetOutput() = res;
  }

  return true;
}

bool DolovVCrsMatMultAll::PostProcessingImpl() {
  return true;
}

}  // namespace dolov_v_crs_mat_mult
