#include "kazennova_a_fox_algorithm/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <thread>
#include <vector>

#include "kazennova_a_fox_algorithm/common/include/common.hpp"

namespace kazennova_a_fox_algorithm {

namespace {

void GetBlock(const std::vector<double> &mat, int rows, int cols, int block_row, int block_col, int block_size,
              double *block_buf) {
  const int start_row = block_row * block_size;
  const int start_col = block_col * block_size;
  const int end_row = std::min(start_row + block_size, rows);
  const int end_col = std::min(start_col + block_size, cols);

  for (int i = 0; i < block_size; ++i) {
    for (int j = 0; j < block_size; ++j) {
      block_buf[(i * block_size) + j] = 0.0;
    }
  }
  for (int i = start_row; i < end_row; ++i) {
    for (int j = start_col; j < end_col; ++j) {
      block_buf[((i - start_row) * block_size) + (j - start_col)] = mat[(i * cols) + j];
    }
  }
}

void MultiplyBlock(const std::vector<double> &block_a, const std::vector<double> &block_b, int block_size, int max_i,
                   int max_j, int max_k, int bi, int bj, int n, std::vector<double> &local_c) {
  for (int i = 0; i < max_i; ++i) {
    const int local_row = (bi * block_size) + i;
    for (int j = 0; j < max_j; ++j) {
      double sum = 0.0;
      for (int kk = 0; kk < max_k; ++kk) {
        sum += block_a[(i * block_size) + kk] * block_b[(kk * block_size) + j];
      }
      local_c[(local_row * n) + ((bj * block_size) + j)] += sum;
    }
  }
}

void RunLocalMultiplication(const std::vector<double> &a, const std::vector<double> &b, int rows_total, int cols_a,
                            int cols_b, int start_row, int local_rows, int block_size, std::vector<double> &local_c) {
  const int blocks_i_local = (local_rows + block_size - 1) / block_size;
  const int blocks_j = (cols_b + block_size - 1) / block_size;
  const int blocks_k = (cols_a + block_size - 1) / block_size;

  int num_threads = static_cast<int>(std::thread::hardware_concurrency());
  if (num_threads <= 0) {
    num_threads = 2;
  }

  std::vector<std::thread> threads;
  threads.reserve(static_cast<size_t>(num_threads));

  std::atomic<size_t> next_block_idx(0);
  const size_t total_blocks = static_cast<size_t>(blocks_i_local) * blocks_j;

  auto worker = [&]() {
    std::vector<double> block_a(static_cast<size_t>(block_size) * block_size);
    std::vector<double> block_b(static_cast<size_t>(block_size) * block_size);

    while (true) {
      const size_t idx = next_block_idx.fetch_add(1);
      if (idx >= total_blocks) {
        break;
      }

      const int bi = static_cast<int>(idx / blocks_j);
      const int bj = static_cast<int>(idx % blocks_j);

      for (int bk = 0; bk < blocks_k; ++bk) {
        const int bi_global = (start_row / block_size) + bi;
        GetBlock(a, rows_total, cols_a, bi_global, bk, block_size, block_a.data());
        GetBlock(b, cols_a, cols_b, bk, bj, block_size, block_b.data());

        const int offset = start_row % block_size;
        const int max_i = std::min(block_size, local_rows - (bi * block_size) - offset);
        const int max_j = std::min(block_size, cols_b - (bj * block_size));
        const int max_k = std::min(block_size, cols_a - (bk * block_size));

        MultiplyBlock(block_a, block_b, block_size, max_i, max_j, max_k, bi, bj, cols_b, local_c);
      }
    }
  };

  for (int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    threads.emplace_back(worker);
  }
  for (auto &thr : threads) {
    thr.join();
  }
}

void ComputeRecvCounts(int world_size, int rows_per_proc, int remainder, int cols_b, std::vector<int> &recv_counts,
                       std::vector<int> &displs, int &total_elements) {
  total_elements = 0;
  for (int proc = 0; proc < world_size; ++proc) {
    const int proc_local_rows = rows_per_proc + (proc < remainder ? 1 : 0);
    recv_counts[proc] = proc_local_rows * cols_b;
    displs[proc] = total_elements;
    total_elements += recv_counts[proc];
  }
}

void GatherAndAssemble(int rank, int world_size, int rows_per_proc, int remainder, int cols_b,
                       const std::vector<double> &local_c, std::vector<double> &c) {
  if (rank == 0) {
    std::vector<int> recv_counts(world_size);
    std::vector<int> displs(world_size);
    int total_elements = 0;
    ComputeRecvCounts(world_size, rows_per_proc, remainder, cols_b, recv_counts, displs, total_elements);
    std::vector<double> gathered(static_cast<size_t>(total_elements));
    MPI_Gatherv(local_c.data(), static_cast<int>(local_c.size()), MPI_DOUBLE, gathered.data(), recv_counts.data(),
                displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for (int proc = 0; proc < world_size; ++proc) {
      const int proc_start_row = (proc * rows_per_proc) + std::min(proc, remainder);
      const int proc_local_rows = rows_per_proc + (proc < remainder ? 1 : 0);
      for (int i = 0; i < proc_local_rows; ++i) {
        for (int j = 0; j < cols_b; ++j) {
          c[((proc_start_row + i) * cols_b) + j] = gathered[displs[proc] + (i * cols_b) + j];
        }
      }
    }
  } else {
    std::vector<int> recv_counts(world_size);
    std::vector<int> displs(world_size);
    int total_elements = 0;
    ComputeRecvCounts(world_size, rows_per_proc, remainder, cols_b, recv_counts, displs, total_elements);
    MPI_Gatherv(local_c.data(), static_cast<int>(local_c.size()), MPI_DOUBLE, nullptr, nullptr, nullptr, MPI_DOUBLE, 0,
                MPI_COMM_WORLD);
  }
}

}  // namespace

KazennovaATestTaskALL::KazennovaATestTaskALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool KazennovaATestTaskALL::ValidationImpl() {
  const auto &in = GetInput();
  if (in.A.data.empty() || in.B.data.empty()) {
    return false;
  }
  if (in.A.rows <= 0 || in.A.cols <= 0 || in.B.rows <= 0 || in.B.cols <= 0) {
    return false;
  }
  if (in.A.cols != in.B.rows) {
    return false;
  }
  return true;
}

bool KazennovaATestTaskALL::PreProcessingImpl() {
  const auto &in = GetInput();
  auto &out = GetOutput();
  out.rows = in.A.rows;
  out.cols = in.B.cols;
  out.data.assign(static_cast<size_t>(out.rows) * out.cols, 0.0);
  return true;
}

bool KazennovaATestTaskALL::RunImpl() {
  int rank = -1;
  int world_size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  const auto &in = GetInput();
  auto &out = GetOutput();

  const int rows_total = in.A.rows;
  const int cols_a = in.A.cols;
  const int cols_b = in.B.cols;
  const auto &a = in.A.data;
  const auto &b = in.B.data;
  auto &c = out.data;

  const int block_size = kBlockSize;

  const int rows_per_proc = rows_total / world_size;
  const int remainder = rows_total % world_size;
  const int start_row = (rank * rows_per_proc) + std::min(rank, remainder);
  const int local_rows = rows_per_proc + (rank < remainder ? 1 : 0);

  if (local_rows == 0) {
    if (rank == 0) {
      return true;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    return true;
  }

  std::vector<double> local_c(static_cast<size_t>(local_rows) * cols_b, 0.0);

  RunLocalMultiplication(a, b, rows_total, cols_a, cols_b, start_row, local_rows, block_size, local_c);

  GatherAndAssemble(rank, world_size, rows_per_proc, remainder, cols_b, local_c, c);

  MPI_Barrier(MPI_COMM_WORLD);
  return true;
}

bool KazennovaATestTaskALL::PostProcessingImpl() {
  return !GetOutput().data.empty();
}

}  // namespace kazennova_a_fox_algorithm
