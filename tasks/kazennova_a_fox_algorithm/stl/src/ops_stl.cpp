#include "kazennova_a_fox_algorithm/stl/include/ops_stl.hpp"

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <thread>
#include <vector>

#include "kazennova_a_fox_algorithm/common/include/common.hpp"
#include "util/include/util.hpp"

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
                   int max_j, int max_k, int bi, int bj, int n, std::vector<double> &c) {
  for (int i = 0; i < max_i; ++i) {
    const int global_row = (bi * block_size) + i;
    for (int j = 0; j < max_j; ++j) {
      const int global_col = (bj * block_size) + j;
      double sum = 0.0;
      for (int kk = 0; kk < max_k; ++kk) {
        sum += block_a[(i * block_size) + kk] * block_b[(kk * block_size) + j];
      }
      c[(global_row * n) + global_col] += sum;
    }
  }
}

}  // namespace

KazennovaATestTaskSTL::KazennovaATestTaskSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool KazennovaATestTaskSTL::ValidationImpl() {
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

bool KazennovaATestTaskSTL::PreProcessingImpl() {
  const auto &in = GetInput();
  auto &out = GetOutput();
  out.rows = in.A.rows;
  out.cols = in.B.cols;
  out.data.assign(static_cast<size_t>(out.rows) * out.cols, 0.0);
  return true;
}

bool KazennovaATestTaskSTL::RunImpl() {
  const auto &in = GetInput();
  auto &out = GetOutput();

  const int m = in.A.rows;
  const int k = in.A.cols;
  const int n = in.B.cols;
  const auto &a = in.A.data;
  const auto &b = in.B.data;
  auto &c = out.data;

  const int bs = kBlockSize;

  const int blocks_i = (m + bs - 1) / bs;
  const int blocks_j = (n + bs - 1) / bs;
  const int blocks_k = (k + bs - 1) / bs;

  int num_threads = ppc::util::GetNumThreads();
  if (num_threads <= 0) {
    num_threads = static_cast<int>(std::thread::hardware_concurrency());
  }
  if (num_threads <= 0) {
    num_threads = 2;
  }

  std::vector<std::thread> threads;
  threads.reserve(static_cast<size_t>(num_threads));

  std::atomic<size_t> next_block_idx(0);
  const size_t total_blocks = static_cast<size_t>(blocks_i) * blocks_j;

  auto worker = [&]() {
    std::vector<double> block_a(static_cast<size_t>(bs) * bs);
    std::vector<double> block_b(static_cast<size_t>(bs) * bs);

    while (true) {
      const size_t idx = next_block_idx.fetch_add(1);
      if (idx >= total_blocks) {
        break;
      }

      const int bi = static_cast<int>(idx / blocks_j);
      const int bj = static_cast<int>(idx % blocks_j);

      for (int bk = 0; bk < blocks_k; ++bk) {
        GetBlock(a, m, k, bi, bk, bs, block_a.data());
        GetBlock(b, k, n, bk, bj, bs, block_b.data());

        const int max_i = std::min(bs, m - (bi * bs));
        const int max_j = std::min(bs, n - (bj * bs));
        const int max_k = std::min(bs, k - (bk * bs));

        MultiplyBlock(block_a, block_b, bs, max_i, max_j, max_k, bi, bj, n, c);
      }
    }
  };

  for (int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    threads.emplace_back(worker);
  }
  for (auto &thr : threads) {
    thr.join();
  }

  return true;
}

bool KazennovaATestTaskSTL::PostProcessingImpl() {
  return !GetOutput().data.empty();
}

}  // namespace kazennova_a_fox_algorithm
