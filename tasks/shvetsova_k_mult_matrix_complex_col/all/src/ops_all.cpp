#include "../../all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <complex>
#include <utility>
#include <vector>

#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/parallel_for.h"
#include "shvetsova_k_mult_matrix_complex_col/common/include/common.hpp"

namespace shvetsova_k_mult_matrix_complex_col {

struct SparseColumn {
  std::vector<int> rows;
  std::vector<std::complex<double>> vals;
};

namespace {

void ComputeColumnTask(int col_idx, const MatrixCCS &matrix_a, const MatrixCCS &matrix_b,
                       std::vector<std::complex<double>> &column_c_local, SparseColumn &out_col) {
  std::ranges::fill(column_c_local, std::complex<double>{0.0, 0.0});

  for (int j = matrix_b.col_ptr[col_idx]; j < matrix_b.col_ptr[col_idx + 1]; j++) {
    int tmp_ind = matrix_b.row_ind[j];
    std::complex<double> tmp_val = matrix_b.values[j];
    for (int ind = matrix_a.col_ptr[tmp_ind]; ind < matrix_a.col_ptr[tmp_ind + 1]; ind++) {
      int row = matrix_a.row_ind[ind];
      std::complex<double> val_a = matrix_a.values[ind];
      column_c_local[row] += tmp_val * val_a;
    }
  }
  for (int index = 0; std::cmp_less(index, column_c_local.size()); ++index) {
    if (column_c_local[index].real() != 0.0 || column_c_local[index].imag() != 0.0) {
      out_col.rows.push_back(index);
      out_col.vals.push_back(column_c_local[index]);
    }
  }
}

void FlattenLocalColumns(const std::vector<SparseColumn> &local_columns, std::vector<int> &local_rows,
                         std::vector<double> &local_vals_real, std::vector<double> &local_vals_imag) {
  for (const auto &col : local_columns) {
    local_rows.insert(local_rows.end(), col.rows.begin(), col.rows.end());
    for (const auto &val : col.vals) {
      local_vals_real.push_back(val.real());
      local_vals_imag.push_back(val.imag());
    }
  }
}

void PrepareGatherCounts(int size, int cols_per_rank, int remainder, int cols, std::vector<int> &recv_counts,
                         std::vector<int> &displs_cols, std::vector<int> &all_nnz_per_col) {
  for (int i = 0; i < size; ++i) {
    recv_counts[i] = cols_per_rank + (i < remainder ? 1 : 0);
    displs_cols[i] = (i * cols_per_rank) + std::min(i, remainder);
  }
  all_nnz_per_col.resize(cols);
}

void PrepareGatherDisplacements(int size, int cols_per_rank, int remainder, const std::vector<int> &all_nnz_per_col,
                                std::vector<int> &recv_counts, std::vector<int> &displs_cols,
                                std::vector<int> &recv_nnz_counts, std::vector<int> &displs_nnz,
                                int &total_global_nnz) {
  for (int i = 0; i < size; ++i) {
    recv_counts[i] = cols_per_rank + (i < remainder ? 1 : 0);
    displs_cols[i] = (i * cols_per_rank) + std::min(i, remainder);

    int nnz_sum = 0;
    for (int j = 0; j < recv_counts[i]; ++j) {
      nnz_sum += all_nnz_per_col[displs_cols[i] + j];
    }
    recv_nnz_counts[i] = nnz_sum;

    if (i > 0) {
      displs_nnz[i] = displs_nnz[i - 1] + recv_nnz_counts[i - 1];
    }
  }
  if (size > 0) {
    total_global_nnz = displs_nnz.back() + recv_nnz_counts.back();
  }
}

// сборка финальной CCS матрицы на 0-м ранге
void AssembleFinalMatrix(const std::vector<int> &all_nnz_per_col, const std::vector<int> &all_rows,
                         const std::vector<double> &all_vals_real, const std::vector<double> &all_vals_imag, int cols,
                         MatrixCCS &matrix_c) {
  matrix_c.row_ind.clear();
  matrix_c.values.clear();
  matrix_c.col_ptr.clear();
  matrix_c.col_ptr.push_back(0);

  int current_nnz = 0;
  for (int i = 0; i < cols; ++i) {
    int nnz = all_nnz_per_col[i];
    for (int k = 0; k < nnz; ++k) {
      matrix_c.row_ind.push_back(all_rows[current_nnz + k]);
      matrix_c.values.emplace_back(all_vals_real[current_nnz + k], all_vals_imag[current_nnz + k]);
    }
    current_nnz += nnz;
    matrix_c.col_ptr.push_back(static_cast<int>(matrix_c.row_ind.size()));
  }
}

}  // namespace

ShvetsovaKMultMatrixComplexALL::ShvetsovaKMultMatrixComplexALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = MatrixCCS(0, 0, std::vector<int>{0}, std::vector<int>{}, std::vector<std::complex<double>>{});
}

bool ShvetsovaKMultMatrixComplexALL::ValidationImpl() {
  return true;
}

bool ShvetsovaKMultMatrixComplexALL::PreProcessingImpl() {
  const auto &matrix_a = std::get<0>(GetInput());
  const auto &matrix_b = std::get<1>(GetInput());

  auto &matrix_c = GetOutput();
  matrix_c.rows = matrix_a.rows;
  matrix_c.cols = matrix_b.cols;
  matrix_c.row_ind.clear();
  matrix_c.values.clear();
  matrix_c.col_ptr.clear();
  matrix_c.col_ptr.push_back(0);
  return true;
}

bool ShvetsovaKMultMatrixComplexALL::RunImpl() {
  const MatrixCCS &matrix_a = std::get<0>(GetInput());
  const MatrixCCS &matrix_b = std::get<1>(GetInput());
  auto &matrix_c = GetOutput();

  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int cols_per_rank = matrix_b.cols / size;
  int remainder = matrix_b.cols % size;

  int local_start = (rank * cols_per_rank) + std::min(rank, remainder);
  int local_count = cols_per_rank + (rank < remainder ? 1 : 0);

  std::vector<SparseColumn> local_columns(local_count);

  tbb::parallel_for(tbb::blocked_range<int>(0, local_count), [&](const tbb::blocked_range<int> &r) {
    std::vector<std::complex<double>> column_c_local(matrix_a.rows, {0.0, 0.0});
    for (int i = r.begin(); i != r.end(); ++i) {
      int global_col = local_start + i;
      ComputeColumnTask(global_col, matrix_a, matrix_b, column_c_local, local_columns[i]);
    }
  });

  std::vector<int> local_nnz_per_col(local_count);
  int total_local_nnz = 0;
  for (int i = 0; i < local_count; ++i) {
    local_nnz_per_col[i] = static_cast<int>(local_columns[i].rows.size());
    total_local_nnz += local_nnz_per_col[i];
  }

  std::vector<int> local_rows;
  local_rows.reserve(total_local_nnz);
  std::vector<double> local_vals_real;
  local_vals_real.reserve(total_local_nnz);
  std::vector<double> local_vals_imag;
  local_vals_imag.reserve(total_local_nnz);

  FlattenLocalColumns(local_columns, local_rows, local_vals_real, local_vals_imag);

  std::vector<int> recv_counts(size, 0);
  std::vector<int> displs_cols(size, 0);
  std::vector<int> all_nnz_per_col;

  if (rank == 0) {
    PrepareGatherCounts(size, cols_per_rank, remainder, matrix_b.cols, recv_counts, displs_cols, all_nnz_per_col);
  }

  MPI_Gatherv(local_nnz_per_col.data(), local_count, MPI_INT, all_nnz_per_col.data(), recv_counts.data(),
              displs_cols.data(), MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<int> recv_nnz_counts(size, 0);
  std::vector<int> displs_nnz(size, 0);
  int total_global_nnz = 0;

  if (rank == 0) {
    PrepareGatherDisplacements(size, cols_per_rank, remainder, all_nnz_per_col, recv_counts, displs_cols,
                               recv_nnz_counts, displs_nnz, total_global_nnz);
  }

  std::vector<int> all_rows;
  std::vector<double> all_vals_real;
  std::vector<double> all_vals_imag;

  if (rank == 0) {
    all_rows.resize(total_global_nnz);
    all_vals_real.resize(total_global_nnz);
    all_vals_imag.resize(total_global_nnz);
  }

  MPI_Gatherv(local_rows.data(), total_local_nnz, MPI_INT, all_rows.data(), recv_nnz_counts.data(), displs_nnz.data(),
              MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Gatherv(local_vals_real.data(), total_local_nnz, MPI_DOUBLE, all_vals_real.data(), recv_nnz_counts.data(),
              displs_nnz.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Gatherv(local_vals_imag.data(), total_local_nnz, MPI_DOUBLE, all_vals_imag.data(), recv_nnz_counts.data(),
              displs_nnz.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    AssembleFinalMatrix(all_nnz_per_col, all_rows, all_vals_real, all_vals_imag, matrix_b.cols, matrix_c);
  }

  return true;
}

bool ShvetsovaKMultMatrixComplexALL::PostProcessingImpl() {
  return true;
}

}  // namespace shvetsova_k_mult_matrix_complex_col
