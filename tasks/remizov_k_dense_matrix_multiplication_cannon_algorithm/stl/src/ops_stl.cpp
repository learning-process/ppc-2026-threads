#include "remizov_k_dense_matrix_multiplication_cannon_algorithm/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstddef>
#include <thread>
#include <utility>
#include <vector>
#include <functional>

#include "remizov_k_dense_matrix_multiplication_cannon_algorithm/common/include/common.hpp"

namespace remizov_k_dense_matrix_multiplication_cannon_algorithm {

template <typename IndexType, typename Func>
static void ParallelFor(IndexType begin, IndexType end, Func&& func) {
  const std::size_t num_threads = std::max(1u, std::thread::hardware_concurrency());
  const IndexType range_length = end - begin;
  if (range_length <= 0) return;

  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  IndexType chunk_size = (range_length + num_threads - 1) / num_threads;
  IndexType start = begin;

  for (std::size_t t = 0; t < num_threads; ++t) {
    IndexType chunk_end = std::min(end, start + chunk_size);
    if (start >= chunk_end) break;

    threads.emplace_back([start, chunk_end, &func]() {
      for (IndexType i = start; i < chunk_end; ++i) {
        func(i);
      }
    });
    start = chunk_end;
  }

  for (auto& th : threads) {
    if (th.joinable()) th.join();
  }
}

template <typename Func>
static void ParallelFor2D(int rows_begin, int rows_end, int cols_begin, int cols_end, Func&& func) {
  const int rows = rows_end - rows_begin;
  const int cols = cols_end - cols_begin;
  const int total = rows * cols;
  if (total <= 0) return;

  ParallelFor(0, total, [&](int linear_idx) {
    int i = rows_begin + linear_idx / cols;
    int j = cols_begin + linear_idx % cols;
    func(i, j);
  });
}


RemizovKDenseMatrixMultiplicationCannonAlgorithmStl::RemizovKDenseMatrixMultiplicationCannonAlgorithmStl(
    const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}


bool RemizovKDenseMatrixMultiplicationCannonAlgorithmStl::ValidationImpl() {
  const auto &input_data = GetInput();
  int block_dim = std::get<0>(input_data);
  const auto &mat_a = std::get<1>(input_data);
  const auto &mat_b = std::get<2>(input_data);

  if (block_dim <= 0) return false;
  if (mat_a.empty() || mat_b.empty()) return false;

  size_t n = mat_a.size();
  if (n != mat_a[0].size()) return false;
  if (n != mat_b.size() || n != mat_b[0].size()) return false;

  return (n % static_cast<size_t>(block_dim) == 0);
}

bool RemizovKDenseMatrixMultiplicationCannonAlgorithmStl::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}


}  // namespace remizov_k_dense_matrix_multiplication_cannon_algorithm
