#include "liulin_y_complex_ccs/all/include/ops_all.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <utility>
#include <vector>

namespace liulin_y_complex_ccs {

namespace {
constexpr double kEpsilon = 1e-15;

std::complex<double> DotProduct(const CCSMatrix &AT, int row_idx, const CCSMatrix &B, int col_idx) {
  std::complex<double> sum(0.0, 0.0);
  int a_ptr = AT.col_index[static_cast<size_t>(row_idx)];
  int a_end = AT.col_index[static_cast<size_t>(row_idx) + 1];
  int b_ptr = B.col_index[static_cast<size_t>(col_idx)];
  int b_end = B.col_index[static_cast<size_t>(col_idx) + 1];

  while (a_ptr < a_end && b_ptr < b_end) {
    int idx_a = AT.row_index[static_cast<size_t>(a_ptr)];
    int idx_b = B.row_index[static_cast<size_t>(b_ptr)];
    if (idx_a == idx_b) {
      sum += AT.values[static_cast<size_t>(a_ptr)] * B.values[static_cast<size_t>(b_ptr)];
      a_ptr++;
      b_ptr++;
    } else if (idx_a < idx_b) {
      a_ptr++;
    } else {
      b_ptr++;
    }
  }
  return sum;
}
}  // namespace

LiulinYComplexCcsAll::LiulinYComplexCcsAll(const InType &in) : Task(in) {
  SetTypeOfTask(GetStaticTypeOfTask());
}

bool LiulinYComplexCcsAll::ValidationImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    const auto &A = GetInput().first;
    const auto &B = GetInput().second;
    if (A.count_cols != B.count_rows) {
      return false;
    }
    if (A.count_rows <= 0 || A.count_cols <= 0 || B.count_cols <= 0) {
      return false;
    }
    if (A.col_index.size() != static_cast<size_t>(A.count_cols + 1)) {
      return false;
    }
    if (B.col_index.size() != static_cast<size_t>(B.count_cols + 1)) {
      return false;
    }
    if (A.col_index[0] != 0 || B.col_index[0] != 0) {
      return false;
    }
  }
  return true;
}

bool LiulinYComplexCcsAll::PreProcessingImpl() {
  return true;
}

void LiulinYComplexCcsAll::BroadcastMatrix(CCSMatrix &mat, int root, MPI_Comm comm) {
  MPI_Bcast(&mat.count_rows, 1, MPI_INT, root, comm);
  MPI_Bcast(&mat.count_cols, 1, MPI_INT, root, comm);
  int nnz = static_cast<int>(mat.values.size());
  MPI_Bcast(&nnz, 1, MPI_INT, root, comm);

  if (mat.col_index.size() != static_cast<size_t>(mat.count_cols + 1)) {
    mat.col_index.resize(static_cast<size_t>(mat.count_cols + 1));
  }
  if (mat.values.size() != static_cast<size_t>(nnz)) {
    mat.values.resize(static_cast<size_t>(nnz));
    mat.row_index.resize(static_cast<size_t>(nnz));
  }

  MPI_Bcast(mat.col_index.data(), static_cast<int>(mat.col_index.size()), MPI_INT, root, comm);
  MPI_Bcast(mat.row_index.data(), static_cast<int>(mat.row_index.size()), MPI_INT, root, comm);

  std::vector<double> real_p(mat.values.size()), imag_p(mat.values.size());
  int rank;
  MPI_Comm_rank(comm, &rank);
  if (rank == root) {
    for (size_t i = 0; i < mat.values.size(); ++i) {
      real_p[i] = mat.values[i].real();
      imag_p[i] = mat.values[i].imag();
    }
  }
  MPI_Bcast(real_p.data(), static_cast<int>(real_p.size()), MPI_DOUBLE, root, comm);
  MPI_Bcast(imag_p.data(), static_cast<int>(imag_p.size()), MPI_DOUBLE, root, comm);
  if (rank != root) {
    for (size_t i = 0; i < mat.values.size(); ++i) {
      mat.values[i] = {real_p[i], imag_p[i]};
    }
  }
}

void LiulinYComplexCcsAll::Transpose(const CCSMatrix &in, CCSMatrix &out) {
  out.count_rows = in.count_cols;
  out.count_cols = in.count_rows;
  out.col_index.assign(static_cast<size_t>(out.count_cols + 1), 0);
  for (int idx : in.row_index) {
    out.col_index[static_cast<size_t>(idx) + 1]++;
  }
  for (size_t i = 0; i < static_cast<size_t>(out.count_cols); ++i) {
    out.col_index[i + 1] += out.col_index[i];
  }

  std::vector<int> current_pos = out.col_index;
  out.row_index.resize(in.values.size());
  out.values.resize(in.values.size());
  for (int j = 0; j < in.count_cols; ++j) {
    for (int k = in.col_index[static_cast<size_t>(j)]; k < in.col_index[static_cast<size_t>(j) + 1]; ++k) {
      int row = in.row_index[static_cast<size_t>(k)];
      int pos = current_pos[static_cast<size_t>(row)]++;
      out.row_index[static_cast<size_t>(pos)] = j;
      out.values[static_cast<size_t>(pos)] = in.values[static_cast<size_t>(k)];
    }
  }
}

bool LiulinYComplexCcsAll::RunImpl() {
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  CCSMatrix A = GetInput().first;
  CCSMatrix B = GetInput().second;

  try {
    BroadcastMatrix(A, 0, MPI_COMM_WORLD);
    BroadcastMatrix(B, 0, MPI_COMM_WORLD);

    CCSMatrix AT;
    Transpose(A, AT);

    std::vector<std::vector<int>> local_rows(static_cast<size_t>(B.count_cols));
    std::vector<std::vector<std::complex<double>>> local_vals(static_cast<size_t>(B.count_cols));
    std::vector<int> local_nnz_per_col(static_cast<size_t>(B.count_cols), 0);

#pragma omp parallel for schedule(dynamic)
    for (int j = rank; j < B.count_cols; j += size) {
      for (int i = 0; i < A.count_rows; ++i) {
        std::complex<double> res = DotProduct(AT, i, B, j);
        if (std::abs(res.real()) > kEpsilon || std::abs(res.imag()) > kEpsilon) {
          local_rows[static_cast<size_t>(j)].push_back(i);
          local_vals[static_cast<size_t>(j)].push_back(res);
          local_nnz_per_col[static_cast<size_t>(j)]++;
        }
      }
    }

    std::vector<int> global_nnz_per_col(static_cast<size_t>(B.count_cols));
    MPI_Allreduce(local_nnz_per_col.data(), global_nnz_per_col.data(), B.count_cols, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    CCSMatrix &C = GetOutput();
    C.count_rows = A.count_rows;
    C.count_cols = B.count_cols;
    C.col_index.assign(static_cast<size_t>(C.count_cols + 1), 0);
    for (size_t i = 0; i < static_cast<size_t>(C.count_cols); ++i) {
      C.col_index[i + 1] = C.col_index[i] + global_nnz_per_col[i];
    }
    int total_nnz = C.col_index.back();
    C.row_index.resize(static_cast<size_t>(total_nnz));
    C.values.resize(static_cast<size_t>(total_nnz));

    for (int j = 0; j < B.count_cols; ++j) {
      int count = global_nnz_per_col[static_cast<size_t>(j)];
      if (count == 0) {
        continue;
      }

      int sender = j % size;
      std::vector<double> r_buf(count), i_buf(count);
      std::vector<int> row_buf(count);

      if (rank == sender) {
        row_buf = local_rows[static_cast<size_t>(j)];
        for (size_t k = 0; k < static_cast<size_t>(count); ++k) {
          r_buf[k] = local_vals[static_cast<size_t>(j)][k].real();
          i_buf[k] = local_vals[static_cast<size_t>(j)][k].imag();
        }
      }

      MPI_Bcast(row_buf.data(), count, MPI_INT, sender, MPI_COMM_WORLD);
      MPI_Bcast(r_buf.data(), count, MPI_DOUBLE, sender, MPI_COMM_WORLD);
      MPI_Bcast(i_buf.data(), count, MPI_DOUBLE, sender, MPI_COMM_WORLD);

      int start_pos = C.col_index[static_cast<size_t>(j)];
      for (int k = 0; k < count; ++k) {
        C.row_index[static_cast<size_t>(start_pos + k)] = row_buf[static_cast<size_t>(k)];
        C.values[static_cast<size_t>(start_pos + k)] = {r_buf[static_cast<size_t>(k)], i_buf[static_cast<size_t>(k)]};
      }
    }
  } catch (...) {
    return false;
  }
  return true;
}

bool LiulinYComplexCcsAll::PostProcessingImpl() {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    const auto &C = GetOutput();
    if (C.count_rows <= 0 || C.count_cols <= 0) {
      return false;
    }
    if (C.col_index.back() != static_cast<int>(C.values.size())) {
      return false;
    }
  }
  return true;
}

}  // namespace liulin_y_complex_ccs
