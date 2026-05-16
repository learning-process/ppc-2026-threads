#include "cheremkhin_a_matr_mult_cannon_alg/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "cheremkhin_a_matr_mult_cannon_alg/common/include/common.hpp"
#include "util/include/util.hpp"

namespace cheremkhin_a_matr_mult_cannon_alg {

namespace {

inline std::size_t Idx(std::size_t n, std::size_t r, std::size_t c) {
  return (r * n) + c;
}

std::size_t CeilDiv(std::size_t a, std::size_t b) {
  return (a + b - 1) / b;
}

int ChooseGridSize(int world_size) {
  if (world_size <= 1) {
    return 1;
  }

  const auto root = static_cast<int>(std::sqrt(static_cast<double>(world_size)));
  for (int grid_dim = root; grid_dim >= 1; --grid_dim) {
    if ((world_size % grid_dim) == 0) {
      return grid_dim;
    }
  }
  return 1;
}

bool IsPerfectSquare(int value) {
  if (value <= 0) {
    return false;
  }
  const int root = static_cast<int>(std::sqrt(static_cast<double>(value)));
  return root * root == value;
}

void MatMulFallbackOMP(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c, std::size_t n) {
  const auto n64 = static_cast<std::int64_t>(n);
#pragma omp parallel for default(none) schedule(static) shared(a, b, c, n, n64)
  for (std::int64_t i = 0; i < n64; ++i) {
    const auto row = static_cast<std::size_t>(i);
    for (std::size_t k = 0; k < n; ++k) {
      const double aik = a[Idx(n, row, k)];
      const double *b_row = b.data() + (k * n);
      double *c_row = c.data() + (row * n);
      for (std::size_t j = 0; j < n; ++j) {
        c_row[j] += aik * b_row[j];
      }
    }
  }
}

void CopyGlobalToPadded(const std::vector<double> &src, std::vector<double> &dst, std::size_t src_n, std::size_t dst_n) {
  const auto src_n64 = static_cast<std::int64_t>(src_n);
#pragma omp parallel for default(none) schedule(static) shared(src, dst, src_n, dst_n, src_n64)
  for (std::int64_t i = 0; i < src_n64; ++i) {
    for (std::size_t j = 0; j < src_n; ++j) {
      dst[Idx(dst_n, static_cast<std::size_t>(i), j)] = src[Idx(src_n, static_cast<std::size_t>(i), j)];
    }
  }
}

void CopyPaddedToGlobal(const std::vector<double> &src, std::vector<double> &dst, std::size_t src_n, std::size_t dst_n) {
  const auto dst_n64 = static_cast<std::int64_t>(dst_n);
#pragma omp parallel for default(none) schedule(static) shared(src, dst, src_n, dst_n, dst_n64)
  for (std::int64_t i = 0; i < dst_n64; ++i) {
    for (std::size_t j = 0; j < dst_n; ++j) {
      dst[Idx(dst_n, static_cast<std::size_t>(i), j)] = src[Idx(src_n, static_cast<std::size_t>(i), j)];
    }
  }
}

void ExtractLocalBlock(const std::vector<double> &src, std::vector<double> &block, std::size_t global_n, std::size_t block_n,
                       int block_row, int block_col) {
  const std::size_t row0 = static_cast<std::size_t>(block_row) * block_n;
  const std::size_t col0 = static_cast<std::size_t>(block_col) * block_n;
  const auto block_n64 = static_cast<std::int64_t>(block_n);
#pragma omp parallel for default(none) schedule(static) shared(src, block, global_n, block_n, row0, col0, block_n64)
  for (std::int64_t i = 0; i < block_n64; ++i) {
    const std::size_t src_row = (row0 + static_cast<std::size_t>(i)) * global_n;
    const std::size_t dst_row = static_cast<std::size_t>(i) * block_n;
    for (std::size_t j = 0; j < block_n; ++j) {
      block[dst_row + j] = src[src_row + col0 + j];
    }
  }
}

void InsertLocalBlock(const std::vector<double> &block, std::vector<double> &dst, std::size_t global_n, std::size_t block_n,
                      int block_row, int block_col) {
  const std::size_t row0 = static_cast<std::size_t>(block_row) * block_n;
  const std::size_t col0 = static_cast<std::size_t>(block_col) * block_n;
  const auto block_n64 = static_cast<std::int64_t>(block_n);
#pragma omp parallel for default(none) schedule(static) shared(block, dst, global_n, block_n, row0, col0, block_n64)
  for (std::int64_t i = 0; i < block_n64; ++i) {
    const std::size_t src_row = static_cast<std::size_t>(i) * block_n;
    const std::size_t dst_row = (row0 + static_cast<std::size_t>(i)) * global_n;
    for (std::size_t j = 0; j < block_n; ++j) {
      dst[dst_row + col0 + j] = block[src_row + j];
    }
  }
}

void ScatterBlocksFromRoot(const std::vector<double> &global_matrix, std::vector<double> &local_block, std::size_t global_n,
                           std::size_t block_n, MPI_Comm cart_comm) {
  int world_rank = 0;
  int world_size = 0;
  MPI_Comm_rank(cart_comm, &world_rank);
  MPI_Comm_size(cart_comm, &world_size);

  if (world_rank == 0) {
    for (int rank = 0; rank < world_size; ++rank) {
      std::array<int, 2> coords = {0, 0};
      MPI_Cart_coords(cart_comm, rank, 2, coords.data());

      if (rank == 0) {
        ExtractLocalBlock(global_matrix, local_block, global_n, block_n, coords[0], coords[1]);
        continue;
      }

      std::vector<double> tmp(block_n * block_n, 0.0);
      ExtractLocalBlock(global_matrix, tmp, global_n, block_n, coords[0], coords[1]);
      MPI_Send(tmp.data(), static_cast<int>(tmp.size()), MPI_DOUBLE, rank, 0, cart_comm);
    }
  } else {
    MPI_Recv(local_block.data(), static_cast<int>(local_block.size()), MPI_DOUBLE, 0, 0, cart_comm, MPI_STATUS_IGNORE);
  }
}

void GatherBlocksToAll(const std::vector<double> &local_block, std::vector<double> &global_matrix, std::size_t global_n,
                       std::size_t block_n, MPI_Comm cart_comm) {
  int world_rank = 0;
  int world_size = 0;
  MPI_Comm_rank(cart_comm, &world_rank);
  MPI_Comm_size(cart_comm, &world_size);

  if (world_rank == 0) {
    for (int rank = 0; rank < world_size; ++rank) {
      std::array<int, 2> coords = {0, 0};
      MPI_Cart_coords(cart_comm, rank, 2, coords.data());

      if (rank == 0) {
        InsertLocalBlock(local_block, global_matrix, global_n, block_n, coords[0], coords[1]);
        continue;
      }

      std::vector<double> tmp(block_n * block_n, 0.0);
      MPI_Recv(tmp.data(), static_cast<int>(tmp.size()), MPI_DOUBLE, rank, 1, cart_comm, MPI_STATUS_IGNORE);
      InsertLocalBlock(tmp, global_matrix, global_n, block_n, coords[0], coords[1]);
    }
  } else {
    MPI_Send(local_block.data(), static_cast<int>(local_block.size()), MPI_DOUBLE, 0, 1, cart_comm);
  }

  MPI_Bcast(global_matrix.data(), static_cast<int>(global_matrix.size()), MPI_DOUBLE, 0, cart_comm);
}

void MulAddLocal(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c, std::size_t block_n) {
  const auto block_n64 = static_cast<std::int64_t>(block_n);

#pragma omp parallel for default(none) schedule(static) shared(a, b, c, block_n, block_n64)
  for (std::int64_t ii = 0; ii < block_n64; ++ii) {
    const auto row = static_cast<std::size_t>(ii);
    const std::size_t a_row = row * block_n;
    const std::size_t c_row = row * block_n;
    double *c_block = c.data() + c_row;
    for (std::size_t kk = 0; kk < block_n; ++kk) {
      const double aik = a[a_row + kk];
      const double *b_block = b.data() + (kk * block_n);
      for (std::int64_t jj = 0; jj < block_n64; ++jj) {
        c_block[jj] += aik * b_block[jj];
      }
    }
  }
}

}  // namespace

CheremkhinAMatrMultCannonAlgALL::CheremkhinAMatrMultCannonAlgALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool CheremkhinAMatrMultCannonAlgALL::ValidationImpl() {
  const std::size_t n = std::get<0>(GetInput());
  const auto &a = std::get<1>(GetInput());
  const auto &b = std::get<2>(GetInput());
  return n > 0 && a.size() == n * n && b.size() == n * n;
}

bool CheremkhinAMatrMultCannonAlgALL::PreProcessingImpl() {
  GetOutput() = {};
  return true;
}

bool CheremkhinAMatrMultCannonAlgALL::RunImpl() {
  const std::size_t n = std::get<0>(GetInput());
  const auto &a_in = std::get<1>(GetInput());
  const auto &b_in = std::get<2>(GetInput());
  const int requested_threads = ppc::util::GetNumThreads();
  int world_rank = 0;
  int world_size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  omp_set_num_threads(requested_threads);

  if (!IsPerfectSquare(world_size)) {
    std::vector<double> out(n * n, 0.0);
    MatMulFallbackOMP(a_in, b_in, out, n);
    GetOutput() = std::move(out);
    return true;
  }

  const int q = ChooseGridSize(world_size);
  const std::size_t block_n = CeilDiv(n, static_cast<std::size_t>(q));
  const std::size_t padded_n = block_n * static_cast<std::size_t>(q);

  std::vector<double> a_padded;
  std::vector<double> b_padded;
  if (world_rank == 0) {
    a_padded.assign(padded_n * padded_n, 0.0);
    b_padded.assign(padded_n * padded_n, 0.0);
    CopyGlobalToPadded(a_in, a_padded, n, padded_n);
    CopyGlobalToPadded(b_in, b_padded, n, padded_n);
  }

  std::vector<double> a_local(block_n * block_n, 0.0);
  std::vector<double> b_local(block_n * block_n, 0.0);
  std::vector<double> c_local(block_n * block_n, 0.0);

  std::array<int, 2> dims = {q, q};
  std::array<int, 2> periods = {1, 1};
  MPI_Comm cart_comm = MPI_COMM_NULL;
  MPI_Cart_create(MPI_COMM_WORLD, 2, dims.data(), periods.data(), 0, &cart_comm);

  std::array<int, 2> coords = {0, 0};
  MPI_Cart_coords(cart_comm, world_rank, 2, coords.data());

  ScatterBlocksFromRoot(a_padded, a_local, padded_n, block_n, cart_comm);
  ScatterBlocksFromRoot(b_padded, b_local, padded_n, block_n, cart_comm);

  int left_rank = 0;
  int right_rank = 0;
  int up_rank = 0;
  int down_rank = 0;
  MPI_Cart_shift(cart_comm, 1, -1, &right_rank, &left_rank);
  MPI_Cart_shift(cart_comm, 0, -1, &down_rank, &up_rank);

  for (int shift = 0; shift < coords[0]; ++shift) {
    MPI_Sendrecv_replace(a_local.data(), static_cast<int>(a_local.size()), MPI_DOUBLE, left_rank, 10, right_rank, 10,
                         cart_comm, MPI_STATUS_IGNORE);
  }
  for (int shift = 0; shift < coords[1]; ++shift) {
    MPI_Sendrecv_replace(b_local.data(), static_cast<int>(b_local.size()), MPI_DOUBLE, up_rank, 11, down_rank, 11,
                         cart_comm, MPI_STATUS_IGNORE);
  }

  for (int step = 0; step < q; ++step) {
    MulAddLocal(a_local, b_local, c_local, block_n);
    if (step + 1 < q) {
      MPI_Sendrecv_replace(a_local.data(), static_cast<int>(a_local.size()), MPI_DOUBLE, left_rank, 20, right_rank, 20,
                           cart_comm, MPI_STATUS_IGNORE);
      MPI_Sendrecv_replace(b_local.data(), static_cast<int>(b_local.size()), MPI_DOUBLE, up_rank, 21, down_rank, 21,
                           cart_comm, MPI_STATUS_IGNORE);
    }
  }

  std::vector<double> c_padded(padded_n * padded_n, 0.0);
  GatherBlocksToAll(c_local, c_padded, padded_n, block_n, cart_comm);

  std::vector<double> out(n * n, 0.0);
  CopyPaddedToGlobal(c_padded, out, padded_n, n);
  MPI_Comm_free(&cart_comm);

  GetOutput() = std::move(out);
  return true;
}

bool CheremkhinAMatrMultCannonAlgALL::PostProcessingImpl() {
  return true;
}

}  // namespace  cheremkhin_a_matr_mult_cannon_alg
