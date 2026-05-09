#include "safronov_m_multiplication_matrix_blockscheme_cannon/stl/include/ops_stl.hpp"

#include <functional>
#include <thread>
#include <utility>
#include <vector>

#include "safronov_m_multiplication_matrix_blockscheme_cannon/common/include/common.hpp"

namespace safronov_m_multiplication_matrix_blocksscheme_cannon {

SafronovMMultiplicationMatrixBlockSchemeCannonSTL::SafronovMMultiplicationMatrixBlockSchemeCannonSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool SafronovMMultiplicationMatrixBlockSchemeCannonSTL::ValidationImpl() {
  const auto &in = GetInput();
  int size_block = std::get<0>(in);
  const auto &matrix_a = std::get<1>(in);
  const auto &matrix_b = std::get<2>(in);
  return (size_block > 0) && (!matrix_a.empty() && !matrix_b.empty()) && (matrix_a.size() == matrix_a[0].size()) &&
         (matrix_b.size() == matrix_b[0].size()) && (matrix_a.size() == matrix_b.size()) &&
         (matrix_a.size() % size_block == 0);
}

bool SafronovMMultiplicationMatrixBlockSchemeCannonSTL::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}

void SafronovMMultiplicationMatrixBlockSchemeCannonSTL::MultiplyingBlocks(std::vector<std::vector<double>> &block_a,
                                                                          std::vector<std::vector<double>> &block_b,
                                                                          std::vector<std::vector<double>> &block_c,
                                                                          int size_block) {
  for (int i = 0; i < size_block; i++) {
    for (int j = 0; j < size_block; j++) {
      for (int k = 0; k < size_block; k++) {
        block_c[i][j] += block_a[i][k] * block_b[k][j];
      }
    }
  }
}

void SafronovMMultiplicationMatrixBlockSchemeCannonSTL::ShiftALeftTask(
    int start, int end, int columns, std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_a) {
  for (int i = start; i < end; i++) {
    std::vector<std::vector<double>> tmp = std::move(matrix_a[i][0]);
    for (int j = 1; j < columns; j++) {
      matrix_a[i][j - 1] = std::move(matrix_a[i][j]);
    }
    matrix_a[i][columns - 1] = std::move(tmp);
  }
}

void SafronovMMultiplicationMatrixBlockSchemeCannonSTL::ShiftBUpTask(
    int start, int end, int columns, std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_b) {
  for (int i = start; i < end; i++) {
    std::vector<std::vector<double>> tmp = std::move(matrix_b[0][i]);
    for (int j = 1; j < columns; j++) {
      matrix_b[j - 1][i] = std::move(matrix_b[j][i]);
    }
    matrix_b[columns - 1][i] = std::move(tmp);
  }
}

void SafronovMMultiplicationMatrixBlockSchemeCannonSTL::CannonTask(
    int start, int end, int columns_blocks, int size_block,
    std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_a,
    std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_b,
    std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_c) {
  for (int j = start; j < end; j++) {
    for (int k = 0; k < columns_blocks; k++) {
      MultiplyingBlocks(matrix_a[j][k], matrix_b[j][k], matrix_c[j][k], size_block);
    }
  }
}

void SafronovMMultiplicationMatrixBlockSchemeCannonSTL::FillingTask(
    int start, int end, int columns_blocks, int size_block,
    const std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocks_c,
    std::vector<std::vector<double>> &matrix_c) {
  for (int i = start; i < end; i++) {
    for (int j = 0; j < columns_blocks; j++) {
      for (int k = 0; k < size_block; k++) {
        for (int col = 0; col < size_block; col++) {
          matrix_c[(i * size_block) + k][(j * size_block) + col] = matrix_blocks_c[i][j][k][col];
        }
      }
    }
  }
}

void SafronovMMultiplicationMatrixBlockSchemeCannonSTL::FillBlocks(
    int start, int end, int columns_blocks, int size_block, const std::vector<std::vector<double>> &matrix_a,
    const std::vector<std::vector<double>> &matrix_b,
    std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocks_a,
    std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocks_b) {
  for (int i = start; i < end; i++) {
    for (int j = 0; j < columns_blocks; j++) {
      int shift = (i + j) % columns_blocks;
      for (int k = 0; k < size_block; k++) {
        for (int col = 0; col < size_block; col++) {
          matrix_blocks_a[i][j][k][col] = matrix_a[(i * size_block) + k][(shift * size_block) + col];
          matrix_blocks_b[i][j][k][col] = matrix_b[(shift * size_block) + k][(j * size_block) + col];
        }
      }
    }
  }
}

void SafronovMMultiplicationMatrixBlockSchemeCannonSTL::ShiftingBlocksMatrixALeft(
    std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocks_a, int columns) {
  unsigned int hardware_threads = std::thread::hardware_concurrency();
  int num_threads = (hardware_threads == 0) ? 1 : static_cast<int>(hardware_threads);
  std::vector<std::thread> threads;
  int base_size = columns / num_threads;
  int rem = columns % num_threads;
  int start = 0;

  for (int thd_id = 0; thd_id < num_threads; thd_id++) {
    int end = start + base_size + (thd_id < rem ? 1 : 0);
    if (start < end) {
      threads.emplace_back(&SafronovMMultiplicationMatrixBlockSchemeCannonSTL::ShiftALeftTask, start, end, columns,
                           std::ref(matrix_blocks_a));
    }
    start = end;
  }
  for (auto &th : threads) {
    th.join();
  }
}

void SafronovMMultiplicationMatrixBlockSchemeCannonSTL::ShiftingBlocksMatrixBUp(
    std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocks_b, int columns) {
  unsigned int hardware_threads = std::thread::hardware_concurrency();
  int num_threads = (hardware_threads == 0) ? 1 : static_cast<int>(hardware_threads);
  std::vector<std::thread> threads;
  int base_size = columns / num_threads;
  int rem = columns % num_threads;
  int start = 0;

  for (int thd_id = 0; thd_id < num_threads; thd_id++) {
    int end = start + base_size + (thd_id < rem ? 1 : 0);
    if (start < end) {
      threads.emplace_back(&SafronovMMultiplicationMatrixBlockSchemeCannonSTL::ShiftBUpTask, start, end, columns,
                           std::ref(matrix_blocks_b));
    }
    start = end;
  }
  for (auto &th : threads) {
    th.join();
  }
}

void SafronovMMultiplicationMatrixBlockSchemeCannonSTL::AlgorithmCannon(
    std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocks_a,
    std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocks_b,
    std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocks_c, int size_block, int columns_blocks) {
  unsigned int hardware_threads = std::thread::hardware_concurrency();
  int num_threads = (hardware_threads == 0) ? 1 : static_cast<int>(hardware_threads);

  for (int i = 0; i < columns_blocks; i++) {
    std::vector<std::thread> threads;
    int base_size = columns_blocks / num_threads;
    int rem = columns_blocks % num_threads;
    int start = 0;

    for (int thd_id = 0; thd_id < num_threads; thd_id++) {
      int end = start + base_size + (thd_id < rem ? 1 : 0);
      if (start < end) {
        threads.emplace_back(&SafronovMMultiplicationMatrixBlockSchemeCannonSTL::CannonTask, start, end, columns_blocks,
                             size_block, std::ref(matrix_blocks_a), std::ref(matrix_blocks_b),
                             std::ref(matrix_blocks_c));
      }
      start = end;
    }
    for (auto &th : threads) {
      th.join();
    }

    if (i < columns_blocks - 1) {
      ShiftingBlocksMatrixALeft(matrix_blocks_a, columns_blocks);
      ShiftingBlocksMatrixBUp(matrix_blocks_b, columns_blocks);
    }
  }
}

void SafronovMMultiplicationMatrixBlockSchemeCannonSTL::FillingResultingMatrix(
    std::vector<std::vector<std::vector<std::vector<double>>>> &matrix_blocks_c,
    std::vector<std::vector<double>> &matrix_c, int size_block, int columns_blocks) {
  unsigned int hardware_threads = std::thread::hardware_concurrency();
  int num_threads = (hardware_threads == 0) ? 1 : static_cast<int>(hardware_threads);
  std::vector<std::thread> threads;
  int base_size = columns_blocks / num_threads;
  int rem = columns_blocks % num_threads;
  int start = 0;

  for (int thd_id = 0; thd_id < num_threads; thd_id++) {
    int end = start + base_size + (thd_id < rem ? 1 : 0);
    if (start < end) {
      threads.emplace_back(&SafronovMMultiplicationMatrixBlockSchemeCannonSTL::FillingTask, start, end, columns_blocks,
                           size_block, std::ref(matrix_blocks_c), std::ref(matrix_c));
    }
    start = end;
  }
  for (auto &th : threads) {
    th.join();
  }
}

bool SafronovMMultiplicationMatrixBlockSchemeCannonSTL::RunImpl() {
  const auto &in = GetInput();
  int size_block = std::get<0>(in);
  const auto &matrix_a = std::get<1>(in);
  const auto &matrix_b = std::get<2>(in);
  int n = static_cast<int>(matrix_a.size());
  int columns_blocks = n / size_block;

  std::vector<std::vector<std::vector<std::vector<double>>>> matrix_blocks_a(
      columns_blocks,
      std::vector<std::vector<std::vector<double>>>(
          columns_blocks, std::vector<std::vector<double>>(size_block, std::vector<double>(size_block))));
  std::vector<std::vector<std::vector<std::vector<double>>>> matrix_blocks_b(
      columns_blocks,
      std::vector<std::vector<std::vector<double>>>(
          columns_blocks, std::vector<std::vector<double>>(size_block, std::vector<double>(size_block))));
  std::vector<std::vector<std::vector<std::vector<double>>>> matrix_blocks_c(
      columns_blocks,
      std::vector<std::vector<std::vector<double>>>(
          columns_blocks, std::vector<std::vector<double>>(size_block, std::vector<double>(size_block, 0.0))));

  unsigned int hardware_threads = std::thread::hardware_concurrency();
  int num_threads = (hardware_threads == 0) ? 1 : static_cast<int>(hardware_threads);
  std::vector<std::thread> threads;
  int base_size = columns_blocks / num_threads;
  int rem = columns_blocks % num_threads;
  int start = 0;

  for (int thd_id = 0; thd_id < num_threads; thd_id++) {
    int end = start + base_size + (thd_id < rem ? 1 : 0);
    if (start < end) {
      threads.emplace_back(&SafronovMMultiplicationMatrixBlockSchemeCannonSTL::FillBlocks, start, end, columns_blocks,
                           size_block, std::ref(matrix_a), std::ref(matrix_b), std::ref(matrix_blocks_a),
                           std::ref(matrix_blocks_b));
    }
    start = end;
  }
  for (auto &th : threads) {
    th.join();
  }

  AlgorithmCannon(matrix_blocks_a, matrix_blocks_b, matrix_blocks_c, size_block, columns_blocks);

  std::vector<std::vector<double>> matrix_c(n, std::vector<double>(n));
  FillingResultingMatrix(matrix_blocks_c, matrix_c, size_block, columns_blocks);
  GetOutput() = std::move(matrix_c);
  return true;
}

bool SafronovMMultiplicationMatrixBlockSchemeCannonSTL::PostProcessingImpl() {
  return true;
}

}  // namespace safronov_m_multiplication_matrix_blocksscheme_cannon
