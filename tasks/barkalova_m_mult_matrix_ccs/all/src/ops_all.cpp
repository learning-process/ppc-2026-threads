#include "barkalova_m_mult_matrix_ccs/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <exception>
#include <utility>
#include <vector>

#include "barkalova_m_mult_matrix_ccs/common/include/common.hpp"

namespace barkalova_m_mult_matrix_ccs {

BarkalovaMMultMatrixCcsALL::BarkalovaMMultMatrixCcsALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    GetInput() = in;
  }

  GetOutput() = CCSMatrix{};
}

bool BarkalovaMMultMatrixCcsALL::ValidationImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank != 0) {
    return true;
  }

  const auto &[A, B] = GetInput();
  if (A.cols != B.rows) {
    return false;
  }
  if (A.rows <= 0 || A.cols <= 0 || B.rows <= 0 || B.cols <= 0) {
    return false;
  }
  if (A.col_ptrs.size() != static_cast<size_t>(A.cols) + 1 || B.col_ptrs.size() != static_cast<size_t>(B.cols) + 1) {
    return false;
  }
  if (A.col_ptrs.empty() || A.col_ptrs[0] != 0 || B.col_ptrs.empty() || B.col_ptrs[0] != 0) {
    return false;
  }
  if (std::cmp_not_equal(A.nnz, A.values.size()) || std::cmp_not_equal(B.nnz, B.values.size())) {
    return false;
  }
  return true;
}

bool BarkalovaMMultMatrixCcsALL::PreProcessingImpl() {
  return true;
}

namespace {
constexpr double kEpsilon = 1e-10;

void TransponirMatr(const CCSMatrix &a, CCSMatrix &at) {
  at.rows = a.cols;
  at.cols = a.rows;
  at.nnz = a.nnz;

  if (a.nnz == 0) {
    at.values.clear();
    at.row_indices.clear();
    at.col_ptrs.assign(at.cols + 1, 0);
    return;
  }

  std::vector<int> row_count(at.cols, 0);
  for (int i = 0; i < a.nnz; i++) {
    row_count[a.row_indices[i]]++;
  }

  at.col_ptrs.resize(at.cols + 1);
  at.col_ptrs[0] = 0;
  for (int i = 0; i < at.cols; i++) {
    at.col_ptrs[i + 1] = at.col_ptrs[i] + row_count[i];
  }

  at.values.resize(a.nnz);
  at.row_indices.resize(a.nnz);

  std::vector<int> current_pos(at.cols, 0);
  for (int col = 0; col < a.cols; col++) {
    for (int i = a.col_ptrs[col]; i < a.col_ptrs[col + 1]; i++) {
      int row = a.row_indices[i];
      Complex val = a.values[i];
      int pos = at.col_ptrs[row] + current_pos[row];
      at.values[pos] = val;
      at.row_indices[pos] = col;
      current_pos[row]++;
    }
  }
}

void BroadcastMatrix(CCSMatrix &matrix) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::array<int, 4> meta = {matrix.rows, matrix.cols, matrix.nnz, 0};
  if (rank == 0) {
    meta[3] = static_cast<int>(matrix.col_ptrs.size());
  }
  MPI_Bcast(meta.data(), 4, MPI_INT, 0, MPI_COMM_WORLD);

  matrix.rows = meta[0];
  matrix.cols = meta[1];
  matrix.nnz = meta[2];
  int col_ptrs_size = meta[3];

  if (matrix.nnz == 0) {
    matrix.values.clear();
    matrix.row_indices.clear();
    matrix.col_ptrs.assign(col_ptrs_size, 0);
    return;
  }

  if (rank != 0) {
    matrix.values.resize(matrix.nnz);
    matrix.row_indices.resize(matrix.nnz);
    matrix.col_ptrs.resize(col_ptrs_size);
  }

  MPI_Bcast(matrix.col_ptrs.data(), col_ptrs_size, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(matrix.row_indices.data(), matrix.nnz, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<double> values_real(matrix.nnz);
  std::vector<double> values_imag(matrix.nnz);

  if (rank == 0) {
    for (int i = 0; i < matrix.nnz; ++i) {
      values_real[i] = matrix.values[i].real();
      values_imag[i] = matrix.values[i].imag();
    }
  }

  MPI_Bcast(values_real.data(), matrix.nnz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(values_imag.data(), matrix.nnz, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    for (int i = 0; i < matrix.nnz; ++i) {
      matrix.values[i] = Complex(values_real[i], values_imag[i]);
    }
  }
}

bool IsNonZero(const Complex &val) {
  return std::abs(val.real()) > kEpsilon || std::abs(val.imag()) > kEpsilon;
}

Complex ComputeScalarProduct(const CCSMatrix &at, const CCSMatrix &b, int row_a, int col_b) {
  Complex sum = Complex(0.0, 0.0);

  int ks = at.col_ptrs[row_a];
  int ls = b.col_ptrs[col_b];
  int kf = at.col_ptrs[row_a + 1];
  int lf = b.col_ptrs[col_b + 1];

  while ((ks < kf) && (ls < lf)) {
    if (at.row_indices[ks] < b.row_indices[ls]) {
      ks++;
    } else if (at.row_indices[ks] > b.row_indices[ls]) {
      ls++;
    } else {
      sum += at.values[ks] * b.values[ls];
      ks++;
      ls++;
    }
  }

  return sum;
}

void ComputeLocalColumns(int start_col, int local_cols, const CCSMatrix &at, const CCSMatrix &b,
                         std::vector<std::vector<int>> &col_rows, std::vector<std::vector<Complex>> &col_vals) {
#pragma omp parallel for schedule(static) default(none) shared(start_col, local_cols, at, b, col_rows, col_vals)
  for (int j_local = 0; j_local < local_cols; ++j_local) {
    int global_col = start_col + j_local;

    std::vector<int> rows;
    std::vector<Complex> vals;
    rows.reserve(100);
    vals.reserve(100);

    for (int i = 0; i < at.cols; i++) {
      Complex sum = ComputeScalarProduct(at, b, i, global_col);
      if (IsNonZero(sum)) {
        rows.push_back(i);
        vals.push_back(sum);
      }
    }

    col_rows[j_local] = std::move(rows);
    col_vals[j_local] = std::move(vals);
  }
}

void LocVectors(int local_cols, const std::vector<std::vector<int>> &col_rows,
                const std::vector<std::vector<Complex>> &col_vals, std::vector<int> &local_row_indices,
                std::vector<Complex> &local_values) {
  for (int j = 0; j < local_cols; ++j) {
    local_row_indices.insert(local_row_indices.end(), col_rows[j].begin(), col_rows[j].end());
    local_values.insert(local_values.end(), col_vals[j].begin(), col_vals[j].end());
  }
}

void GatherResults(int rank, int size, int local_nnz, const std::vector<int> &local_row_indices,
                   const std::vector<Complex> &local_values, std::vector<int> &global_row_indices,
                   std::vector<double> &global_values_real, std::vector<double> &global_values_imag, int total_nnz) {
  std::vector<int> recv_counts(size, 0);
  MPI_Gather(&local_nnz, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<int> displs(size, 0);
  if (rank == 0) {
    for (int i = 1; i < size; ++i) {
      displs[i] = displs[i - 1] + recv_counts[i - 1];
    }
  }

  if (rank == 0) {
    global_row_indices.resize(total_nnz);
  }
  MPI_Gatherv(local_row_indices.data(), local_nnz, MPI_INT, global_row_indices.data(), recv_counts.data(),
              displs.data(), MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<double> local_values_real(local_nnz);
  std::vector<double> local_values_imag(local_nnz);
  for (int i = 0; i < local_nnz; ++i) {
    local_values_real[i] = local_values[i].real();
    local_values_imag[i] = local_values[i].imag();
  }

  if (rank == 0) {
    global_values_real.resize(total_nnz);
  }
  MPI_Gatherv(local_values_real.data(), local_nnz, MPI_DOUBLE, global_values_real.data(), recv_counts.data(),
              displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    global_values_imag.resize(total_nnz);
  }
  MPI_Gatherv(local_values_imag.data(), local_nnz, MPI_DOUBLE, global_values_imag.data(), recv_counts.data(),
              displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void BroadcastRes(int rank, int total_rows, int total_cols, int total_nnz, std::vector<int> &global_col_ptrs,
                  std::vector<int> &global_row_indices, std::vector<double> &global_values_real,
                  std::vector<double> &global_values_imag) {
  int bcast_rows = total_rows;
  int bcast_cols = total_cols;
  int bcast_nnz = total_nnz;

  MPI_Bcast(&bcast_rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&bcast_cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&bcast_nnz, 1, MPI_INT, 0, MPI_COMM_WORLD);

  int col_ptrs_size = static_cast<int>(global_col_ptrs.size());
  MPI_Bcast(&col_ptrs_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (rank != 0) {
    global_col_ptrs.resize(col_ptrs_size);
  }
  MPI_Bcast(global_col_ptrs.data(), col_ptrs_size, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    global_row_indices.resize(bcast_nnz);
  }
  MPI_Bcast(global_row_indices.data(), bcast_nnz, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    global_values_real.resize(bcast_nnz);
    global_values_imag.resize(bcast_nnz);
  }
  MPI_Bcast(global_values_real.data(), bcast_nnz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(global_values_imag.data(), bcast_nnz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

}  // namespace

bool BarkalovaMMultMatrixCcsALL::RunImpl() {
  try {
    int rank = 0;
    int size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    CCSMatrix a = GetInput().first;
    CCSMatrix b = GetInput().second;

    BroadcastMatrix(a);
    BroadcastMatrix(b);

    CCSMatrix at;
    TransponirMatr(a, at);

    const int total_cols = b.cols;
    const int total_rows = a.rows;

    int cols_per_process = total_cols / size;
    int remainder = total_cols % size;

    int start_col = (rank * cols_per_process) + std::min(rank, remainder);
    int local_cols = cols_per_process + (rank < remainder ? 1 : 0);

    std::vector<std::vector<int>> col_rows(local_cols);
    std::vector<std::vector<Complex>> col_vals(local_cols);

    ComputeLocalColumns(start_col, local_cols, at, b, col_rows, col_vals);

    std::vector<int> local_row_indices;
    std::vector<Complex> local_values;
    LocVectors(local_cols, col_rows, col_vals, local_row_indices, local_values);
    int local_nnz = static_cast<int>(local_values.size());

    std::vector<int> global_col_nnz(total_cols, 0);
    for (int j = 0; j < local_cols; ++j) {
      global_col_nnz[start_col + j] = static_cast<int>(col_rows[j].size());
    }

    std::vector<int> recv_col_nnz(total_cols);
    MPI_Allreduce(global_col_nnz.data(), recv_col_nnz.data(), total_cols, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    std::vector<int> global_col_ptrs = {0};
    for (int j = 0; j < total_cols; ++j) {
      global_col_ptrs.push_back(global_col_ptrs.back() + recv_col_nnz[j]);
    }
    int total_nnz = global_col_ptrs.back();

    std::vector<int> global_row_indices;
    std::vector<double> global_values_real;
    std::vector<double> global_values_imag;
    GatherResults(rank, size, local_nnz, local_row_indices, local_values, global_row_indices, global_values_real,
                  global_values_imag, total_nnz);

    BroadcastRes(rank, total_rows, total_cols, total_nnz, global_col_ptrs, global_row_indices, global_values_real,
                 global_values_imag);

    std::vector<Complex> global_values(total_nnz);
    for (int i = 0; i < total_nnz; ++i) {
      global_values[i] = Complex(global_values_real[i], global_values_imag[i]);
    }

    CCSMatrix c;
    c.rows = total_rows;
    c.cols = total_cols;
    c.nnz = total_nnz;
    c.values = std::move(global_values);
    c.row_indices = std::move(global_row_indices);
    c.col_ptrs = std::move(global_col_ptrs);

    GetOutput() = c;
    return true;

  } catch (const std::exception &) {
    return false;
  }
}

bool BarkalovaMMultMatrixCcsALL::PostProcessingImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank != 0) {
    return true;
  }

  const auto &c = GetOutput();
  if (c.rows <= 0 || c.cols <= 0) {
    return false;
  }
  if (c.col_ptrs.size() != static_cast<size_t>(c.cols) + 1) {
    return false;
  }
  for (size_t i = 1; i < c.col_ptrs.size(); ++i) {
    if (c.col_ptrs[i] < c.col_ptrs[i - 1]) {
      return false;
    }
  }
  if (std::cmp_not_equal(c.nnz, c.values.size()) || std::cmp_not_equal(c.nnz, c.row_indices.size())) {
    return false;
  }
  if (c.col_ptrs[0] != 0) {
    return false;
  }
  if (c.col_ptrs.back() != c.nnz) {
    return false;
  }
  return true;
}

}  // namespace barkalova_m_mult_matrix_ccs
