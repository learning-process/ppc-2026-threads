#include "goriacheva_k_mult_sparse_complex_matrix_ccs/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <ranges>
#include <utility>
#include <vector>

#include "goriacheva_k_mult_sparse_complex_matrix_ccs/common/include/common.hpp"

namespace goriacheva_k_mult_sparse_complex_matrix_ccs {

GoriachevaKMultSparseComplexMatrixCcsALL::GoriachevaKMultSparseComplexMatrixCcsALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool GoriachevaKMultSparseComplexMatrixCcsALL::ValidationImpl() {
  auto &[a, b] = GetInput();
  return (a.cols == b.rows) && !a.col_ptr.empty() && !b.col_ptr.empty();
}

bool GoriachevaKMultSparseComplexMatrixCcsALL::PreProcessingImpl() {
  GetOutput() = {};
  return true;
}

namespace {

void ProcessColumn(int j, const SparseMatrixCCS &a, const SparseMatrixCCS &b, std::vector<Complex> &values,
                   std::vector<int> &rows) {
  std::vector<Complex> accumulator(a.rows);
  std::vector<int> marker(a.rows, -1);
  std::vector<int> used_rows;

  for (int bi = b.col_ptr[j]; bi < b.col_ptr[j + 1]; bi++) {
    int k = b.row_ind[bi];
    Complex b_val = b.values[bi];

    for (int ai = a.col_ptr[k]; ai < a.col_ptr[k + 1]; ai++) {
      int i = a.row_ind[ai];

      if (marker[i] != j) {
        marker[i] = j;
        accumulator[i] = Complex(0.0, 0.0);
        used_rows.push_back(i);
      }

      accumulator[i] += a.values[ai] * b_val;
    }
  }

  std::ranges::sort(used_rows);

  for (int r : used_rows) {
    if (accumulator[r] != Complex(0.0, 0.0)) {
      rows.push_back(r);
      values.push_back(accumulator[r]);
    }
  }
}

void ComputeDispls(const std::vector<int> &recv_counts, std::vector<int> &displs, int &total_nnz) {
  total_nnz = 0;
  for (std::size_t i = 0; i < recv_counts.size(); i++) {
    displs[i] = total_nnz;
    total_nnz += recv_counts[i];
  }
}

void BuildGlobalColPtr(int size, int cols_per_proc, int total_cols, const std::vector<int> &local_col_ptr,
                       std::vector<int> &global_col_ptr, int rank, MPI_Comm comm) {
  int offset = 0;
  int global_col = 0;

  for (int proc = 0; proc < size; proc++) {
    int start = proc * cols_per_proc;
    int end = std::min(start + cols_per_proc, total_cols);

    int p_cols = 0;
    if (start < end) {
      p_cols = end - start;
    }

    std::vector<int> tmp_col_ptr(p_cols + 1);

    if (rank == proc) {
      tmp_col_ptr = local_col_ptr;
    }

    MPI_Bcast(tmp_col_ptr.data(), p_cols + 1, MPI_INT, proc, comm);

    if (rank == 0) {
      for (int j = 0; j < p_cols; j++) {
        global_col_ptr[global_col++] = offset + tmp_col_ptr[j];
      }
      offset += tmp_col_ptr[p_cols];
    }
  }
}

}  // namespace

bool GoriachevaKMultSparseComplexMatrixCcsALL::RunImpl() {
  auto &a = std::get<0>(GetInput());
  auto &b = std::get<1>(GetInput());
  auto &c = GetOutput();

  int rank = 0;
  int size = 1;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  c.rows = a.rows;
  c.cols = b.cols;

  int cols_per_proc = (b.cols + size - 1) / size;
  int start = rank * cols_per_proc;
  int end = std::min(start + cols_per_proc, b.cols);
  int local_cols = 0;
  if (start < end) {
    local_cols = end - start;
  }

  std::vector<std::vector<Complex>> local_values(local_cols);
  std::vector<std::vector<int>> local_rows(local_cols);

#pragma omp parallel for default(none) shared(a, b, local_values, local_rows, start, local_cols)
  for (int j = 0; j < local_cols; j++) {
    ProcessColumn(start + j, a, b, local_values[j], local_rows[j]);
  }

  std::vector<int> local_col_ptr(local_cols + 1);
  int local_nnz = 0;

  for (int j = 0; j < local_cols; j++) {
    local_col_ptr[j] = local_nnz;
    local_nnz += static_cast<int>(local_values[j].size());
  }
  local_col_ptr[local_cols] = local_nnz;

  std::vector<Complex> local_vals;
  std::vector<int> local_inds;

  local_vals.reserve(local_nnz);
  local_inds.reserve(local_nnz);

  for (int j = 0; j < local_cols; j++) {
    local_vals.insert(local_vals.end(), local_values[j].begin(), local_values[j].end());
    local_inds.insert(local_inds.end(), local_rows[j].begin(), local_rows[j].end());
  }

  std::vector<int> recv_counts(size);
  MPI_Gather(&local_nnz, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<int> displs(size, 0);
  int total_nnz = 0;

  if (rank == 0) {
    ComputeDispls(recv_counts, displs, total_nnz);
  }

  std::vector<int> recv_counts_dbl(size);
  std::vector<int> displs_dbl(size);

  if (rank == 0) {
    for (int i = 0; i < size; i++) {
      recv_counts_dbl[i] = recv_counts[i] * 2;
      displs_dbl[i] = displs[i] * 2;
    }
  }

  std::vector<Complex> global_vals;
  std::vector<int> global_inds;

  if (rank == 0) {
    global_vals.resize(total_nnz);
    global_inds.resize(total_nnz);
  }

  MPI_Gatherv(local_vals.data(), local_nnz * 2, MPI_DOUBLE, global_vals.data(), recv_counts_dbl.data(),
              displs_dbl.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Gatherv(local_inds.data(), local_nnz, MPI_INT, global_inds.data(), recv_counts.data(), displs.data(), MPI_INT, 0,
              MPI_COMM_WORLD);

  if (rank == 0) {
    c.col_ptr.resize(c.cols + 1);
  }

  BuildGlobalColPtr(size, cols_per_proc, b.cols, local_col_ptr, c.col_ptr, rank, MPI_COMM_WORLD);

  if (rank == 0) {
    c.col_ptr[c.cols] = total_nnz;
    c.values = std::move(global_vals);
    c.row_ind = std::move(global_inds);
  }

  MPI_Bcast(&total_nnz, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    c.rows = a.rows;
    c.cols = b.cols;
    c.values.resize(total_nnz);
    c.row_ind.resize(total_nnz);
    c.col_ptr.resize(c.cols + 1);
  }

  MPI_Bcast(c.values.data(), total_nnz * 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(c.row_ind.data(), total_nnz, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(c.col_ptr.data(), c.cols + 1, MPI_INT, 0, MPI_COMM_WORLD);

  return true;
}

bool GoriachevaKMultSparseComplexMatrixCcsALL::PostProcessingImpl() {
  return true;
}

}  // namespace goriacheva_k_mult_sparse_complex_matrix_ccs
