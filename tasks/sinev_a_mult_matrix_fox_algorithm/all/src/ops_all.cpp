#include "sinev_a_mult_matrix_fox_algorithm/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include "sinev_a_mult_matrix_fox_algorithm/common/include/common.hpp"

namespace sinev_a_mult_matrix_fox_algorithm {

SinevAMultMatrixFoxAlgorithmALL::SinevAMultMatrixFoxAlgorithmALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());

  GetInput() = in;
  GetOutput() = {};
}

bool SinevAMultMatrixFoxAlgorithmALL::ValidationImpl() {
  const auto &[n, a, b] = GetInput();

  return (n > 0U) && (a.size() == (n * n)) && (b.size() == (n * n));
}

bool SinevAMultMatrixFoxAlgorithmALL::PreProcessingImpl() {
  const auto &[n, a, b] = GetInput();

  GetOutput().resize(n * n, 0.0);

  return true;
}

void SinevAMultMatrixFoxAlgorithmALL::SimpleMultiply(size_t n, const std::vector<double> &a,
                                                     const std::vector<double> &b, std::vector<double> &c) {
#pragma omp parallel for default(none) shared(n, a, b, c) collapse(2)
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      double sum = 0.0;

      for (size_t k = 0; k < n; ++k) {
        sum += a[(i * n) + k] * b[(k * n) + j];
      }

      c[(i * n) + j] = sum;
    }
  }
}

void SinevAMultMatrixFoxAlgorithmALL::DecomposeToBlocks(const std::vector<double> &src, std::vector<double> &dst,
                                                        size_t n, size_t bs, int q) {
#pragma omp parallel for default(none) shared(src, dst, n, bs, q) collapse(2)
  for (int bi = 0; bi < q; ++bi) {
    for (int bj = 0; bj < q; ++bj) {
      const size_t block_offset = static_cast<size_t>((bi * q) + bj) * (bs * bs);

      for (size_t i = 0; i < bs; ++i) {
        for (size_t j = 0; j < bs; ++j) {
          const size_t src_idx = ((static_cast<size_t>(bi) * bs + i) * n) + (static_cast<size_t>(bj) * bs + j);

          const size_t dst_idx = block_offset + (i * bs) + j;

          dst[dst_idx] = src[src_idx];
        }
      }
    }
  }
}

void SinevAMultMatrixFoxAlgorithmALL::AssembleFromBlocks(const std::vector<double> &src, std::vector<double> &dst,
                                                         size_t n, size_t bs, int q) {
#pragma omp parallel for default(none) shared(src, dst, n, bs, q) collapse(2)
  for (int bi = 0; bi < q; ++bi) {
    for (int bj = 0; bj < q; ++bj) {
      const size_t block_offset = static_cast<size_t>((bi * q) + bj) * (bs * bs);

      for (size_t i = 0; i < bs; ++i) {
        for (size_t j = 0; j < bs; ++j) {
          const size_t src_idx = block_offset + (i * bs) + j;

          const size_t dst_idx = ((static_cast<size_t>(bi) * bs + i) * n) + (static_cast<size_t>(bj) * bs + j);

          dst[dst_idx] = src[src_idx];
        }
      }
    }
  }
}

void SinevAMultMatrixFoxAlgorithmALL::LocalMatrixMultiply(const std::vector<double> &local_a,
                                                          const std::vector<double> &local_b,
                                                          std::vector<double> &local_c, size_t bs) {
#pragma omp parallel for default(none) shared(local_a, local_b, local_c, bs) collapse(2)
  for (size_t i = 0; i < bs; ++i) {
    for (size_t j = 0; j < bs; ++j) {
      double sum = 0.0;

      for (size_t k = 0; k < bs; ++k) {
        sum += local_a[(i * bs) + k] * local_b[(k * bs) + j];
      }

      local_c[(i * bs) + j] += sum;
    }
  }
}

bool SinevAMultMatrixFoxAlgorithmALL::NeedFallback(size_t n, int q, int world_size) {
  return ((q * q) != world_size) || ((n % static_cast<size_t>(q)) != 0U);
}

void SinevAMultMatrixFoxAlgorithmALL::ExecuteFallback(int rank, size_t n, const std::vector<double> &a,
                                                      const std::vector<double> &b, std::vector<double> &c) {
  if (rank == 0) {
    SimpleMultiply(n, a, b, c);
  }

  MPI_Bcast(c.data(), static_cast<int>(n * n), MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void SinevAMultMatrixFoxAlgorithmALL::ScatterBlocks(int rank, const std::vector<double> &blocks_a,
                                                    const std::vector<double> &blocks_b, std::vector<double> &local_a,
                                                    std::vector<double> &local_b, size_t block_size) {
  const double *send_a = (rank == 0) ? blocks_a.data() : nullptr;

  const double *send_b = (rank == 0) ? blocks_b.data() : nullptr;

  MPI_Scatter(send_a, static_cast<int>(block_size), MPI_DOUBLE, local_a.data(), static_cast<int>(block_size),
              MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Scatter(send_b, static_cast<int>(block_size), MPI_DOUBLE, local_b.data(), static_cast<int>(block_size),
              MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void SinevAMultMatrixFoxAlgorithmALL::RunFoxStages(int q, int row, int col, size_t bs, size_t block_size,
                                                   MPI_Comm row_comm, std::vector<double> &local_a,
                                                   std::vector<double> &local_b, std::vector<double> &local_c) {
  std::vector<double> temp_a(block_size);

  for (int step = 0; step < q; ++step) {
    const int root = (row + step) % q;

    if (col == root) {
      temp_a = local_a;
    }

    MPI_Bcast(temp_a.data(), static_cast<int>(block_size), MPI_DOUBLE, root, row_comm);

    LocalMatrixMultiply(temp_a, local_b, local_c, bs);

    const int send_to = (((row - 1 + q) % q) * q) + col;

    const int recv_from = (((row + 1) % q) * q) + col;

    MPI_Sendrecv_replace(local_b.data(), static_cast<int>(block_size), MPI_DOUBLE, send_to, 0, recv_from, 0,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
}

void SinevAMultMatrixFoxAlgorithmALL::GatherResult(int rank, int world_size, size_t n, size_t bs, size_t block_size,
                                                   int q, const std::vector<double> &local_c, std::vector<double> &c) {
  std::vector<double> blocks_c;

  if (rank == 0) {
    blocks_c.resize(static_cast<size_t>(world_size) * block_size);
  }

  double *recv_buffer = (rank == 0) ? blocks_c.data() : nullptr;

  MPI_Gather(local_c.data(), static_cast<int>(block_size), MPI_DOUBLE, recv_buffer, static_cast<int>(block_size),
             MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    AssembleFromBlocks(blocks_c, c, n, bs, q);
  }

  MPI_Bcast(c.data(), static_cast<int>(n * n), MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

bool SinevAMultMatrixFoxAlgorithmALL::RunImpl() {
  int rank = 0;
  int world_size = 1;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  const auto &[n, a, b] = GetInput();

  auto &c = GetOutput();

  const int q = static_cast<int>(std::sqrt(world_size));

  if (NeedFallback(n, q, world_size)) {
    ExecuteFallback(rank, n, a, b, c);

    return true;
  }

  const size_t bs = n / static_cast<size_t>(q);

  const size_t block_size = bs * bs;

  const int row = rank / q;
  const int col = rank % q;

  std::vector<double> local_a(block_size);
  std::vector<double> local_b(block_size);
  std::vector<double> local_c(block_size, 0.0);

  std::vector<double> blocks_a;
  std::vector<double> blocks_b;

  if (rank == 0) {
    blocks_a.resize(static_cast<size_t>(world_size) * block_size);

    blocks_b.resize(static_cast<size_t>(world_size) * block_size);

    DecomposeToBlocks(a, blocks_a, n, bs, q);

    DecomposeToBlocks(b, blocks_b, n, bs, q);
  }

  ScatterBlocks(rank, blocks_a, blocks_b, local_a, local_b, block_size);

  MPI_Comm row_comm = MPI_COMM_NULL;

  const int color = row;
  const int key = col;

  MPI_Comm_split(MPI_COMM_WORLD, color, key, &row_comm);

  RunFoxStages(q, row, col, bs, block_size, row_comm, local_a, local_b, local_c);

  GatherResult(rank, world_size, n, bs, block_size, q, local_c, c);

  if (row_comm != MPI_COMM_NULL) {
    MPI_Comm_free(&row_comm);
  }

  return true;
}

bool SinevAMultMatrixFoxAlgorithmALL::PostProcessingImpl() {
  return true;
}

}  // namespace sinev_a_mult_matrix_fox_algorithm
