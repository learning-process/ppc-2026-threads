#include "sinev_a_mult_matrix_fox_algorithm/stl/include/ops_stl.hpp"

#include <atomic>
#include <cmath>
#include <cstddef>
#include <thread>
#include <vector>

#include "sinev_a_mult_matrix_fox_algorithm/common/include/common.hpp"

namespace sinev_a_mult_matrix_fox_algorithm {

SinevAMultMatrixFoxAlgorithmSTL::SinevAMultMatrixFoxAlgorithmSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool SinevAMultMatrixFoxAlgorithmSTL::ValidationImpl() {
  const auto &[matrix_size, matrix_a, matrix_b] = GetInput();
  return matrix_size > 0 && matrix_a.size() == matrix_size * matrix_size &&
         matrix_b.size() == matrix_size * matrix_size;
}

bool SinevAMultMatrixFoxAlgorithmSTL::PreProcessingImpl() {
  const auto &[matrix_size, matrix_a, matrix_b] = GetInput();
  GetOutput() = std::vector<double>(matrix_size * matrix_size, 0.0);
  return true;
}

void SinevAMultMatrixFoxAlgorithmSTL::SimpleMultiply(size_t n, const std::vector<double> &a,
                                                     const std::vector<double> &b, std::vector<double> &c) {
  for (size_t i = 0; i < n; ++i) {
    for (size_t k = 0; k < n; ++k) {
      double tmp = a[(i * n) + k];
      for (size_t j = 0; j < n; ++j) {
        c[(i * n) + j] += tmp * b[(k * n) + j];
      }
    }
  }
}

void SinevAMultMatrixFoxAlgorithmSTL::DecomposeToBlocks(const std::vector<double> &src, std::vector<double> &dst,
                                                        size_t n, size_t bs, int q) {
  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 2;
  }

  std::vector<std::thread> threads;
  threads.reserve(num_threads);
  std::atomic<size_t> next_block(0);
  size_t total_blocks = static_cast<size_t>(q) * static_cast<size_t>(q);

  for (unsigned int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    threads.emplace_back([&]() {
      size_t block_idx = 0;
      while ((block_idx = next_block.fetch_add(1)) < total_blocks) {
        int bi = static_cast<int>(block_idx / q);
        int bj = static_cast<int>(block_idx % q);

        const size_t block_off = block_idx * (bs * bs);
        for (size_t i = 0; i < bs; ++i) {
          for (size_t j = 0; j < bs; ++j) {
            const size_t src_idx = ((static_cast<size_t>(bi) * bs + i) * n) + (static_cast<size_t>(bj) * bs + j);
            const size_t dst_idx = block_off + (i * bs) + j;
            dst[dst_idx] = src[src_idx];
          }
        }
      }
    });
  }

  for (auto &thread : threads) {
    thread.join();
  }
}

void SinevAMultMatrixFoxAlgorithmSTL::AssembleFromBlocks(const std::vector<double> &src, std::vector<double> &dst,
                                                         size_t n, size_t bs, int q) {
  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 2;
  }

  std::vector<std::thread> threads;
  threads.reserve(num_threads);
  std::atomic<size_t> next_block(0);
  size_t total_blocks = static_cast<size_t>(q) * static_cast<size_t>(q);

  for (unsigned int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    threads.emplace_back([&]() {
      size_t block_idx = 0;
      while ((block_idx = next_block.fetch_add(1)) < total_blocks) {
        int bi = static_cast<int>(block_idx / q);
        int bj = static_cast<int>(block_idx % q);

        const size_t block_off = block_idx * (bs * bs);
        for (size_t i = 0; i < bs; ++i) {
          for (size_t j = 0; j < bs; ++j) {
            const size_t src_idx = block_off + (i * bs) + j;
            const size_t dst_idx = ((static_cast<size_t>(bi) * bs + i) * n) + (static_cast<size_t>(bj) * bs + j);
            dst[dst_idx] = src[src_idx];
          }
        }
      }
    });
  }

  for (auto &thread : threads) {
    thread.join();
  }
}

void SinevAMultMatrixFoxAlgorithmSTL::MultiplyBlocks(const std::vector<double> &blocks_a,
                                                     const std::vector<double> &blocks_b, std::vector<double> &blocks_c,
                                                     size_t bs, size_t a_off, size_t b_off, size_t c_off) {
  for (size_t ii = 0; ii < bs; ++ii) {
    for (size_t kk = 0; kk < bs; ++kk) {
      const double val = blocks_a[a_off + (ii * bs) + kk];
      const size_t b_base = b_off + (kk * bs);
      const size_t c_base = c_off + (ii * bs);
      for (size_t jj = 0; jj < bs; ++jj) {
        blocks_c[c_base + jj] += val * blocks_b[b_base + jj];
      }
    }
  }
}

void SinevAMultMatrixFoxAlgorithmSTL::FoxStep(const std::vector<double> &blocks_a, const std::vector<double> &blocks_b,
                                              std::vector<double> &blocks_c, size_t bs, int q, int step) {
  const size_t block_size = bs * bs;
  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 2;
  }

  std::vector<std::thread> threads;
  threads.reserve(num_threads);
  std::atomic<size_t> next_cell(0);
  size_t total_cells = static_cast<size_t>(q) * static_cast<size_t>(q);

  for (unsigned int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    threads.emplace_back([&]() {
      size_t cell_idx = 0;
      while ((cell_idx = next_cell.fetch_add(1)) < total_cells) {
        int i = static_cast<int>(cell_idx / q);
        int j = static_cast<int>(cell_idx % q);
        const int k = (i + step) % q;

        const size_t a_off = (static_cast<size_t>((i * q) + k)) * block_size;
        const size_t b_off = (static_cast<size_t>((k * q) + j)) * block_size;
        const size_t c_off = (static_cast<size_t>((i * q) + j)) * block_size;

        MultiplyBlocks(blocks_a, blocks_b, blocks_c, bs, a_off, b_off, c_off);
      }
    });
  }

  for (auto &thread : threads) {
    thread.join();
  }
}

bool SinevAMultMatrixFoxAlgorithmSTL::RunImpl() {
  const auto &input = GetInput();
  const size_t n = std::get<0>(input);
  const auto &a = std::get<1>(input);
  const auto &b = std::get<2>(input);
  auto &c = GetOutput();

  if (n <= 64) {
    SimpleMultiply(n, a, b, c);
    return true;
  }

  size_t bs = 64;
  while (n % bs != 0 && bs > 16) {
    bs /= 2;
  }

  if (n % bs != 0) {
    SimpleMultiply(n, a, b, c);
    return true;
  }

  const int actual_q = static_cast<int>(n / bs);

  const auto total_blocks = static_cast<size_t>(actual_q) * static_cast<size_t>(actual_q);
  const auto block_elements = bs * bs;

  std::vector<double> blocks_a(total_blocks * block_elements);
  std::vector<double> blocks_b(total_blocks * block_elements);
  std::vector<double> blocks_c(total_blocks * block_elements, 0.0);

  DecomposeToBlocks(a, blocks_a, n, bs, actual_q);
  DecomposeToBlocks(b, blocks_b, n, bs, actual_q);

  for (int step = 0; step < actual_q; ++step) {
    FoxStep(blocks_a, blocks_b, blocks_c, bs, actual_q, step);
  }

  AssembleFromBlocks(blocks_c, c, n, bs, actual_q);

  return true;
}

size_t SinevAMultMatrixFoxAlgorithmSTL::ChooseBlockSize(size_t n) {
  if (n % 128 == 0) {
    return 128;
  }
  if (n % 64 == 0) {
    return 64;
  }
  if (n % 32 == 0) {
    return 32;
  }
  if (n % 16 == 0) {
    return 16;
  }
  return 1;
}

bool SinevAMultMatrixFoxAlgorithmSTL::PostProcessingImpl() {
  return true;
}

}  // namespace sinev_a_mult_matrix_fox_algorithm
