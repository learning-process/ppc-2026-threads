#include "chyokotov_a_dense_matrix_mul_foxs_algorithm/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <thread>
#include <utility>
#include <vector>

#include "chyokotov_a_dense_matrix_mul_foxs_algorithm/common/include/common.hpp"
#include "util/include/util.hpp"

namespace chyokotov_a_dense_matrix_mul_foxs_algorithm {

ChyokotovADenseMatMulFoxAlgorithmSTL::ChyokotovADenseMatMulFoxAlgorithmSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool ChyokotovADenseMatMulFoxAlgorithmSTL::ValidationImpl() {
  return (GetInput().first.size() == GetInput().second.size());
}

bool ChyokotovADenseMatMulFoxAlgorithmSTL::PreProcessingImpl() {
  GetOutput().clear();
  GetOutput().resize(GetInput().first.size(), 0.0);
  return true;
}

int ChyokotovADenseMatMulFoxAlgorithmSTL::CalculateBlockSize(int n) {
  return static_cast<int>(std::sqrt(static_cast<double>(n)));
}

int ChyokotovADenseMatMulFoxAlgorithmSTL::CountBlock(int n, int size) {
  return (n + size - 1) / size;
}

std::vector<std::pair<int, int>> ChyokotovADenseMatMulFoxAlgorithmSTL::Blocks(int count_block) {
  std::vector<std::pair<int, int>> blocks;
  for (int ic = 0; ic < count_block; ic++) {
    for (int jc = 0; jc < count_block; jc++) {
      blocks.emplace_back(ic, jc);
    }
  }
  return blocks;
}

void ChyokotovADenseMatMulFoxAlgorithmSTL::Matmul(std::vector<double> &a, std::vector<double> &b, int n, int istart,
                                                  int iend, int jstart, int jend, int kstart, int kend) {
  for (int i = istart; i < iend; i++) {
    for (int j = jstart; j < jend; j++) {
      double sum = 0.0;
      for (int k = kstart; k < kend; k++) {
        sum += a[(i * n) + k] * b[(k * n) + j];
      }
      GetOutput()[(i * n) + j] += sum;
    }
  }
}

bool ChyokotovADenseMatMulFoxAlgorithmSTL::RunImpl() {
  std::vector<double> a = GetInput().first;
  std::vector<double> b = GetInput().second;
  int n = static_cast<int>(std::sqrt(static_cast<double>(a.size())));
  if (n == 0) {
    return true;
  }

  int block_size = CalculateBlockSize(n);
  int count_block = CountBlock(n, block_size);
  std::vector<std::pair<int, int>> blocks = Blocks(count_block);

  auto num_threads = static_cast<size_t>(ppc::util::GetNumThreads());
  if (num_threads == 0) {
    num_threads = 4;
  }

  std::vector<std::thread> threads;
  size_t blocks_per_thread = blocks.size() / num_threads;

  for (size_t tt = 0; tt < num_threads; ++tt) {
    size_t start_idx = tt * blocks_per_thread;
    size_t end_idx = (tt == num_threads - 1) ? blocks.size() : start_idx + blocks_per_thread;

    threads.emplace_back([&, start_idx, end_idx]() {
      for (size_t idx = start_idx; idx < end_idx; ++idx) {
        auto [ic, jc] = blocks[idx];

        int istart = ic * block_size;
        int jstart = jc * block_size;
        int iend = std::min(istart + block_size, n);
        int jend = std::min(jstart + block_size, n);

        for (int kc = 0; kc < count_block; kc++) {
          int kstart = kc * block_size;
          int kend = std::min(kstart + block_size, n);
          Matmul(a, b, n, istart, iend, jstart, jend, kstart, kend);
        }
      }
    });
  }

  for (auto &thread : threads) {
    thread.join();
  }

  return true;
}

bool ChyokotovADenseMatMulFoxAlgorithmSTL::PostProcessingImpl() {
  return true;
}

}  // namespace chyokotov_a_dense_matrix_mul_foxs_algorithm
