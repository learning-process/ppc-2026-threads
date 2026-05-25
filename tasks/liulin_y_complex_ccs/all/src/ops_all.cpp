#include "liulin_y_complex_ccs/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <cmath>
#include <complex>
#include <cstddef>
#include <utility>
#include <vector>

#include "liulin_y_complex_ccs/common/include/common.hpp"

namespace liulin_y_complex_ccs {

namespace {
constexpr double kEpsilon = 1e-15;

std::complex<double> DotProduct(const CCSMatrix &mat_at, int row_idx, const CCSMatrix &mat_b, int col_idx) {
  std::complex<double> sum(0.0, 0.0);
  int a_ptr = mat_at.col_index[static_cast<size_t>(row_idx)];
  int a_end = mat_at.col_index[static_cast<size_t>(row_idx) + 1];
  int b_ptr = mat_b.col_index[static_cast<size_t>(col_idx)];
  int b_end = mat_b.col_index[static_cast<size_t>(col_idx) + 1];

  while (a_ptr < a_end && b_ptr < b_end) {
    int idx_a = mat_at.row_index[static_cast<size_t>(a_ptr)];
    int idx_b = mat_b.row_index[static_cast<size_t>(b_ptr)];
    if (idx_a == idx_b) {
      sum += mat_at.values[static_cast<size_t>(a_ptr)] * mat_b.values[static_cast<size_t>(b_ptr)];
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

void GatherColumnData(int col_idx, int sender, int count, int rank, const std::vector<int> &local_rows,
                      const std::vector<std::complex<double>> &local_vals, CCSMatrix &mat_c) {
  std::vector<double> r_buf(static_cast<size_t>(count));
  std::vector<double> i_buf(static_cast<size_t>(count));
  std::vector<int> row_buf(static_cast<size_t>(count));

  if (rank == sender) {
    row_buf = local_rows;
    for (size_t k = 0; std::cmp_less(k, count); ++k) {
      r_buf[k] = local_vals[k].real();
      i_buf[k] = local_vals[k].imag();
    }
  }

  MPI_Bcast(row_buf.data(), count, MPI_INT, sender, MPI_COMM_WORLD);
  MPI_Bcast(r_buf.data(), count, MPI_DOUBLE, sender, MPI_COMM_WORLD);
  MPI_Bcast(i_buf.data(), count, MPI_DOUBLE, sender, MPI_COMM_WORLD);

  auto start_pos = static_cast<size_t>(mat_c.col_index[static_cast<size_t>(col_idx)]);
  for (size_t k = 0; std::cmp_less(k, count); ++k) {
    mat_c.row_index[start_pos + k] = row_buf[k];
    mat_c.values[start_pos + k] = {r_buf[k], i_buf[k]};
  }
}

void BroadcastMatrix(CCSMatrix &mat, int root, MPI_Comm comm) {
  MPI_Bcast(&mat.count_rows, 1, MPI_INT, root, comm);
  MPI_Bcast(&mat.count_cols, 1, MPI_INT, root, comm);
  int nnz = static_cast<int>(mat.values.size());
  MPI_Bcast(&nnz, 1, MPI_INT, root, comm);

  if (std::cmp_not_equal(mat.col_index.size(), mat.count_cols + 1)) {
    mat.col_index.resize(static_cast<size_t>(mat.count_cols) + 1);
  }
  if (std::cmp_not_equal(mat.values.size(), nnz)) {
    mat.values.resize(static_cast<size_t>(nnz));
    mat.row_index.resize(static_cast<size_t>(nnz));
  }

  MPI_Bcast(mat.col_index.data(), static_cast<int>(mat.col_index.size()), MPI_INT, root, comm);
  MPI_Bcast(mat.row_index.data(), static_cast<int>(mat.row_index.size()), MPI_INT, root, comm);

  std::vector<double> real_p(mat.values.size());
  std::vector<double> imag_p(mat.values.size());
  int rank = 0;
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

void Transpose(const CCSMatrix &in, CCSMatrix &out) {
  out.count_rows = in.count_cols;
  out.count_cols = in.count_rows;
  out.col_index.assign(static_cast<size_t>(out.count_cols) + 1, 0);
  for (int idx_r : in.row_index) {
    out.col_index[static_cast<size_t>(idx_r) + 1]++;
  }
  for (size_t i = 0; std::cmp_less(i, out.count_cols); ++i) {
    out.col_index[i + 1] += out.col_index[i];
  }

  std::vector<int> current_pos = out.col_index;
  out.row_index.resize(in.values.size());
  out.values.resize(in.values.size());
  for (int idx_j = 0; idx_j < in.count_cols; ++idx_j) {
    for (int idx_k = in.col_index[static_cast<size_t>(idx_j)]; idx_k < in.col_index[static_cast<size_t>(idx_j) + 1];
         ++idx_k) {
      int row = in.row_index[static_cast<size_t>(idx_k)];
      int pos = current_pos[static_cast<size_t>(row)]++;
      out.row_index[static_cast<size_t>(pos)] = idx_j;
      out.values[static_cast<size_t>(pos)] = in.values[static_cast<size_t>(idx_k)];
    }
  }
}
}  // namespace

LiulinYComplexCcsAll::LiulinYComplexCcsAll(const InType &in) : BaseTask() {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool LiulinYComplexCcsAll::ValidationImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    const auto &mat_a = GetInput().first;
    const auto &mat_b = GetInput().second;
    if (mat_a.count_cols != mat_b.count_rows) {
      return false;
    }
    if (mat_a.count_rows <= 0 || mat_a.count_cols <= 0 || mat_b.count_cols <= 0) {
      return false;
    }
    if (std::cmp_not_equal(mat_a.col_index.size(), mat_a.count_cols + 1)) {
      return false;
    }
    if (std::cmp_not_equal(mat_b.col_index.size(), mat_b.count_cols + 1)) {
      return false;
    }
  }
  return true;
}

bool LiulinYComplexCcsAll::PreProcessingImpl() {
  return true;
}

bool LiulinYComplexCcsAll::RunImpl() {
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  CCSMatrix mat_a = GetInput().first;
  CCSMatrix mat_b = GetInput().second;

  try {
    BroadcastMatrix(mat_a, 0, MPI_COMM_WORLD);
    BroadcastMatrix(mat_b, 0, MPI_COMM_WORLD);
    CCSMatrix mat_at;
    Transpose(mat_a, mat_at);

    std::vector<std::vector<int>> local_rows(static_cast<size_t>(mat_b.count_cols));
    std::vector<std::vector<std::complex<double>>> local_vals(static_cast<size_t>(mat_b.count_cols));
    std::vector<int> local_nnz(static_cast<size_t>(mat_b.count_cols), 0);

#pragma omp parallel for default(none) shared(rank, size, mat_a, mat_b, mat_at, local_rows, local_vals, local_nnz) \
    schedule(dynamic)
    for (int idx_j = rank; idx_j < mat_b.count_cols; idx_j += size) {
      for (int idx_i = 0; idx_i < mat_a.count_rows; ++idx_i) {
        std::complex<double> res = DotProduct(mat_at, idx_i, mat_b, idx_j);
        if (std::abs(res.real()) > kEpsilon || std::abs(res.imag()) > kEpsilon) {
          local_rows[static_cast<size_t>(idx_j)].push_back(idx_i);
          local_vals[static_cast<size_t>(idx_j)].push_back(res);
          local_nnz[static_cast<size_t>(idx_j)]++;
        }
      }
    }

    std::vector<int> global_nnz(static_cast<size_t>(mat_b.count_cols));
    MPI_Allreduce(local_nnz.data(), global_nnz.data(), mat_b.count_cols, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    CCSMatrix &mat_res = GetOutput();
    mat_res.count_rows = mat_a.count_rows;
    mat_res.count_cols = mat_b.count_cols;
    mat_res.col_index.assign(static_cast<size_t>(mat_res.count_cols) + 1, 0);
    for (size_t i = 0; std::cmp_less(i, mat_res.count_cols); ++i) {
      mat_res.col_index[i + 1] = mat_res.col_index[i] + global_nnz[i];
    }
    mat_res.row_index.resize(static_cast<size_t>(mat_res.col_index.back()));
    mat_res.values.resize(static_cast<size_t>(mat_res.col_index.back()));

    for (int idx_j = 0; idx_j < mat_b.count_cols; ++idx_j) {
      int count = global_nnz[static_cast<size_t>(idx_j)];
      if (count > 0) {
        GatherColumnData(idx_j, idx_j % size, count, rank, local_rows[static_cast<size_t>(idx_j)],
                         local_vals[static_cast<size_t>(idx_j)], mat_res);
      }
    }
  } catch (...) {
    return false;
  }
  return true;
}

bool LiulinYComplexCcsAll::PostProcessingImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    const auto &mat_res = GetOutput();
    if (mat_res.count_rows <= 0 || mat_res.count_cols <= 0) {
      return false;
    }
  }
  return true;
}

}  // namespace liulin_y_complex_ccs
