#include "luzan_e_double_sparse_matrix_mult/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <cmath>
#include <utility>
#include <vector>

#include "luzan_e_double_sparse_matrix_mult/common/include/common.hpp"

namespace luzan_e_double_sparse_matrix_mult {

LuzanEDoubleSparseMatrixMultALL::LuzanEDoubleSparseMatrixMultALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());

  GetInput() = in;
  // GetOutput() = 0;
}

void LuzanEDoubleSparseMatrixMultALL::BroadcastMatrix(SparseMatrix &m, int root) {
  MPI_Bcast(&m.rows, 1, MPI_UNSIGNED, root, MPI_COMM_WORLD);
  MPI_Bcast(&m.cols, 1, MPI_UNSIGNED, root, MPI_COMM_WORLD);

  int nnz = static_cast<int>(m.value.size());
  int ci_size = static_cast<int>(m.col_index.size());
  MPI_Bcast(&nnz, 1, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Bcast(&ci_size, 1, MPI_INT, root, MPI_COMM_WORLD);

  m.value.resize(nnz);
  m.row.resize(nnz);
  m.col_index.resize(ci_size);

  MPI_Bcast(m.value.data(), nnz, MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Bcast(m.row.data(), nnz, MPI_UNSIGNED, root, MPI_COMM_WORLD);
  MPI_Bcast(m.col_index.data(), ci_size, MPI_UNSIGNED, root, MPI_COMM_WORLD);
}

void LuzanEDoubleSparseMatrixMultALL::BuildColDistribution(int b_cols, int nprocs, std::vector<int> &counts,
                                                           std::vector<int> &displs) {
  counts.resize(nprocs);
  displs.resize(nprocs, 0);

  int base = b_cols / nprocs;
  int rem = b_cols % nprocs;
  for (int i = 0; i < nprocs; i++) {
    counts[i] = base + (i < rem ? 1 : 0);
  }
  for (int i = 1; i < nprocs; i++) {
    displs[i] = displs[i - 1] + counts[i - 1];
  }
}

void LuzanEDoubleSparseMatrixMultALL::ComputeLocalCols(const SparseMatrix &a, const SparseMatrix &b, int col_start,
                                                       int col_count, std::vector<std::vector<double>> &values_per_col,
                                                       std::vector<std::vector<unsigned>> &rows_per_col) {
  values_per_col.resize(col_count);
  rows_per_col.resize(col_count);

#pragma omp parallel for schedule(static) default(none) \
    shared(a, b, values_per_col, rows_per_col, col_start, col_count, kEPS)
  for (int lc = 0; lc < col_count; lc++) {
    int b_col = col_start + lc;

    std::vector<double> tmp_col(a.rows, 0.0);

    unsigned b_rs = b.col_index[b_col];
    unsigned b_re = b.col_index[b_col + 1];

    for (unsigned bp = b_rs; bp < b_re; bp++) {
      double b_val = b.value[bp];
      unsigned b_row = b.row[bp];

      unsigned a_rs = a.col_index[b_row];
      unsigned a_re = a.col_index[b_row + 1];

      for (unsigned ap = a_rs; ap < a_re; ap++) {
        tmp_col[a.row[ap]] += a.value[ap] * b_val;
      }
    }

    for (unsigned i = 0; i < a.rows; i++) {
      if (std::fabs(tmp_col[i]) > kEPS) {
        values_per_col[lc].push_back(tmp_col[i]);
        rows_per_col[lc].push_back(i);
      }
    }
  }
}

void LuzanEDoubleSparseMatrixMultALL::FlattenLocalCols(const std::vector<std::vector<double>> &values_per_col,
                                                       const std::vector<std::vector<unsigned>> &rows_per_col,
                                                       std::vector<int> &col_nnz, std::vector<double> &flat_vals,
                                                       std::vector<unsigned> &flat_rows) {
  int col_count = static_cast<int>(values_per_col.size());
  col_nnz.resize(col_count);

  for (int lc = 0; lc < col_count; lc++) {
    col_nnz[lc] = static_cast<int>(values_per_col[lc].size());
    flat_vals.insert(flat_vals.end(), values_per_col[lc].begin(), values_per_col[lc].end());
    flat_rows.insert(flat_rows.end(), rows_per_col[lc].begin(), rows_per_col[lc].end());
  }
}

void LuzanEDoubleSparseMatrixMultALL::GatherFlatArrays(int rank, int nprocs, const std::vector<double> &local_vals,
                                                       const std::vector<unsigned> &local_rows,
                                                       std::vector<double> &global_vals,
                                                       std::vector<unsigned> &global_rows) {
  int local_nnz = static_cast<int>(local_vals.size());

  std::vector<int> nnz_counts(nprocs, 0);
  std::vector<int> nnz_displs(nprocs, 0);
  MPI_Gather(&local_nnz, 1, MPI_INT, nnz_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    for (int i = 1; i < nprocs; i++) {
      nnz_displs[i] = nnz_displs[i - 1] + nnz_counts[i - 1];
    }
    int total_nnz = nnz_displs[nprocs - 1] + nnz_counts[nprocs - 1];
    global_vals.resize(total_nnz);
    global_rows.resize(total_nnz);
  }

  MPI_Gatherv(local_vals.data(), local_nnz, MPI_DOUBLE, global_vals.data(), nnz_counts.data(), nnz_displs.data(),
              MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Gatherv(local_rows.data(), local_nnz, MPI_UNSIGNED, global_rows.data(), nnz_counts.data(), nnz_displs.data(),
              MPI_UNSIGNED, 0, MPI_COMM_WORLD);
}

SparseMatrix LuzanEDoubleSparseMatrixMultALL::CalcProdMPIOMP(const SparseMatrix &a_in, const SparseMatrix &b_in) {
  int rank = 0;
  int nprocs = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  SparseMatrix a = (rank == 0) ? a_in : SparseMatrix{};
  SparseMatrix b = (rank == 0) ? b_in : SparseMatrix{};
  BroadcastMatrix(a, 0);
  BroadcastMatrix(b, 0);

  int b_cols = static_cast<int>(b.cols);
  std::vector<int> counts;
  std::vector<int> displs;
  BuildColDistribution(b_cols, nprocs, counts, displs);

  int my_col_start = displs[rank];
  int my_col_count = counts[rank];

  std::vector<std::vector<double>> values_per_col;
  std::vector<std::vector<unsigned>> rows_per_col;
  ComputeLocalCols(a, b, my_col_start, my_col_count, values_per_col, rows_per_col);

  std::vector<int> col_nnz;
  std::vector<double> flat_vals;
  std::vector<unsigned> flat_rows;
  FlattenLocalCols(values_per_col, rows_per_col, col_nnz, flat_vals, flat_rows);

  std::vector<int> global_col_nnz(rank == 0 ? b_cols : 1);
  MPI_Gatherv(col_nnz.data(), my_col_count, MPI_INT, global_col_nnz.data(), counts.data(), displs.data(), MPI_INT, 0,
              MPI_COMM_WORLD);

  std::vector<double> global_vals;
  std::vector<unsigned> global_rows;
  GatherFlatArrays(rank, nprocs, flat_vals, flat_rows, global_vals, global_rows);

  SparseMatrix c;
  if (rank == 0) {
    c.rows = a.rows;
    c.cols = b.cols;
    c.value = std::move(global_vals);
    c.row = std::move(global_rows);

    c.col_index.reserve(static_cast<unsigned>(b_cols) + 1);
    c.col_index.push_back(0);
    for (int j = 0; j < b_cols; j++) {
      c.col_index.push_back(c.col_index.back() + static_cast<unsigned>(global_col_nnz[j]));
    }
  }

  return c;
}

bool LuzanEDoubleSparseMatrixMultALL::ValidationImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    return true;
  }
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());
  return a.GetCols() == b.GetRows() && a.GetCols() != 0 && a.GetRows() != 0 && b.GetCols() != 0;
}

bool LuzanEDoubleSparseMatrixMultALL::PreProcessingImpl() {
  return true;
}

bool LuzanEDoubleSparseMatrixMultALL::RunImpl() {
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());

  GetOutput() = CalcProdMPIOMP(a, b);
  return true;
}

bool LuzanEDoubleSparseMatrixMultALL::PostProcessingImpl() {
  return true;
}

}  // namespace luzan_e_double_sparse_matrix_mult
