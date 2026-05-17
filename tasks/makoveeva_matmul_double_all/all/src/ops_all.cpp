#include "makoveeva_matmul_double_all/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace makoveeva_matmul_double_all {

MatmulDoubleAllTask::MatmulDoubleAllTask(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<double>();
}

bool MatmulDoubleAllTask::ValidationImpl() {
  const auto &input = GetInput();
  const size_t n = std::get<0>(input);
  const auto &a = std::get<1>(input);
  const auto &b = std::get<2>(input);

  return n > 0 && a.size() == n * n && b.size() == n * n;
}

bool MatmulDoubleAllTask::PreProcessingImpl() {
  const auto &input = GetInput();
  matrix_size_ = std::get<0>(input);
  matrix_a_ = std::get<1>(input);
  matrix_b_ = std::get<2>(input);
  result_matrix_.assign(matrix_size_ * matrix_size_, 0.0);

  return true;
}

void MatmulDoubleAllTask::ParallelMultiply(size_t n, const std::vector<double> &a, const std::vector<double> &b,
                                           std::vector<double> &c) {
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

void MatmulDoubleAllTask::SplitIntoBlocks(const std::vector<double> &src, std::vector<double> &dst, size_t n, size_t bs,
                                          int grid_size) {
#pragma omp parallel for default(none) shared(src, dst, n, bs, grid_size) collapse(2)
  for (int bi = 0; bi < grid_size; ++bi) {
    for (int bj = 0; bj < grid_size; ++bj) {
      const size_t block_start = static_cast<size_t>((bi * grid_size) + bj) * (bs * bs);

      for (size_t i = 0; i < bs; ++i) {
        for (size_t j = 0; j < bs; ++j) {
          const size_t src_pos = ((static_cast<size_t>(bi) * bs + i) * n) + (static_cast<size_t>(bj) * bs + j);
          const size_t dst_pos = block_start + (i * bs) + j;
          dst[dst_pos] = src[src_pos];
        }
      }
    }
  }
}

void MatmulDoubleAllTask::MergeFromBlocks(const std::vector<double> &src, std::vector<double> &dst, size_t n, size_t bs,
                                          int grid_size) {
#pragma omp parallel for default(none) shared(src, dst, n, bs, grid_size) collapse(2)
  for (int bi = 0; bi < grid_size; ++bi) {
    for (int bj = 0; bj < grid_size; ++bj) {
      const size_t block_start = static_cast<size_t>((bi * grid_size) + bj) * (bs * bs);

      for (size_t i = 0; i < bs; ++i) {
        for (size_t j = 0; j < bs; ++j) {
          const size_t src_pos = block_start + (i * bs) + j;
          const size_t dst_pos = ((static_cast<size_t>(bi) * bs + i) * n) + (static_cast<size_t>(bj) * bs + j);
          dst[dst_pos] = src[src_pos];
        }
      }
    }
  }
}

void MatmulDoubleAllTask::MultiplyBlockPair(const std::vector<double> &block_a, const std::vector<double> &block_b,
                                            std::vector<double> &block_c, size_t bs) {
#pragma omp parallel for default(none) shared(block_a, block_b, block_c, bs) collapse(2)
  for (size_t i = 0; i < bs; ++i) {
    for (size_t j = 0; j < bs; ++j) {
      double sum = 0.0;
      for (size_t k = 0; k < bs; ++k) {
        sum += block_a[(i * bs) + k] * block_b[(k * bs) + j];
      }
      block_c[(i * bs) + j] += sum;
    }
  }
}

bool MatmulDoubleAllTask::IsValidConfiguration(size_t n, int grid_size, int world_size) {
  return ((grid_size * grid_size) == world_size) && ((n % static_cast<size_t>(grid_size)) == 0);
}

void MatmulDoubleAllTask::HandleFallback(int rank, size_t n, const std::vector<double> &a, const std::vector<double> &b,
                                         std::vector<double> &c) {
  if (rank == 0) {
    ParallelMultiply(n, a, b, c);
  }
  MPI_Bcast(c.data(), static_cast<int>(n * n), MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void MatmulDoubleAllTask::DistributeBlocks(int rank, const std::vector<double> &blocks_a,
                                           const std::vector<double> &blocks_b, std::vector<double> &local_a,
                                           std::vector<double> &local_b, size_t block_sz) {
  const double *send_a = (rank == 0) ? blocks_a.data() : nullptr;
  const double *send_b = (rank == 0) ? blocks_b.data() : nullptr;

  MPI_Scatter(send_a, static_cast<int>(block_sz), MPI_DOUBLE, local_a.data(), static_cast<int>(block_sz), MPI_DOUBLE, 0,
              MPI_COMM_WORLD);

  MPI_Scatter(send_b, static_cast<int>(block_sz), MPI_DOUBLE, local_b.data(), static_cast<int>(block_sz), MPI_DOUBLE, 0,
              MPI_COMM_WORLD);
}

void MatmulDoubleAllTask::ExecuteFoxIterations(int grid, int row_id, int col_id, size_t bs, size_t block_sz,
                                               MPI_Comm row_comm, std::vector<double> &local_a,
                                               std::vector<double> &local_b, std::vector<double> &local_c) {
  std::vector<double> broadcast_buffer(block_sz);

  for (int stage = 0; stage < grid; ++stage) {
    const int source = (row_id + stage) % grid;

    if (col_id == source) {
      broadcast_buffer = local_a;
    }

    MPI_Bcast(broadcast_buffer.data(), static_cast<int>(block_sz), MPI_DOUBLE, source, row_comm);

    MultiplyBlockPair(broadcast_buffer, local_b, local_c, bs);

    const int target = (((row_id - 1 + grid) % grid) * grid) + col_id;
    const int origin = (((row_id + 1) % grid) * grid) + col_id;

    MPI_Sendrecv_replace(local_b.data(), static_cast<int>(block_sz), MPI_DOUBLE, target, 0, origin, 0, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);
  }
}

void MatmulDoubleAllTask::CollectResults(int rank, int world_size, size_t n, size_t bs, size_t block_sz, int grid,
                                         const std::vector<double> &local_c, std::vector<double> &c) {
  std::vector<double> all_blocks;

  if (rank == 0) {
    all_blocks.resize(static_cast<size_t>(world_size) * block_sz);
  }

  double *recv_buf = (rank == 0) ? all_blocks.data() : nullptr;

  MPI_Gather(local_c.data(), static_cast<int>(block_sz), MPI_DOUBLE, recv_buf, static_cast<int>(block_sz), MPI_DOUBLE,
             0, MPI_COMM_WORLD);

  if (rank == 0) {
    MergeFromBlocks(all_blocks, c, n, bs, grid);
  }

  MPI_Bcast(c.data(), static_cast<int>(n * n), MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

bool MatmulDoubleAllTask::RunImpl() {
  int process_rank = 0;
  int total_processes = 1;

  MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &total_processes);

  const size_t n = matrix_size_;
  const auto &a = matrix_a_;
  const auto &b = matrix_b_;
  auto &c = result_matrix_;

  const int grid_size = static_cast<int>(std::sqrt(total_processes));

  if (!IsValidConfiguration(n, grid_size, total_processes)) {
    HandleFallback(process_rank, n, a, b, c);
    GetOutput() = c;
    return true;
  }

  const size_t bs = n / static_cast<size_t>(grid_size);
  const size_t block_sz = bs * bs;

  const int row_idx = process_rank / grid_size;
  const int col_idx = process_rank % grid_size;

  std::vector<double> local_a_block(block_sz);
  std::vector<double> local_b_block(block_sz);
  std::vector<double> local_c_block(block_sz, 0.0);

  std::vector<double> all_blocks_a;
  std::vector<double> all_blocks_b;

  if (process_rank == 0) {
    all_blocks_a.resize(static_cast<size_t>(total_processes) * block_sz);
    all_blocks_b.resize(static_cast<size_t>(total_processes) * block_sz);

    SplitIntoBlocks(a, all_blocks_a, n, bs, grid_size);
    SplitIntoBlocks(b, all_blocks_b, n, bs, grid_size);
  }

  DistributeBlocks(process_rank, all_blocks_a, all_blocks_b, local_a_block, local_b_block, block_sz);

  MPI_Comm row_communicator = MPI_COMM_NULL;
  MPI_Comm_split(MPI_COMM_WORLD, row_idx, col_idx, &row_communicator);

  ExecuteFoxIterations(grid_size, row_idx, col_idx, bs, block_sz, row_communicator, local_a_block, local_b_block,
                       local_c_block);

  CollectResults(process_rank, total_processes, n, bs, block_sz, grid_size, local_c_block, c);

  if (row_communicator != MPI_COMM_NULL) {
    MPI_Comm_free(&row_communicator);
  }

  GetOutput() = c;
  return true;
}

bool MatmulDoubleAllTask::PostProcessingImpl() {
  return true;
}

}  // namespace makoveeva_matmul_double_all
