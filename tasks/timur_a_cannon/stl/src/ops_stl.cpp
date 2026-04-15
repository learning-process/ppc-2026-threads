#include "timur_a_cannon/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstddef>
#include <thread>
#include <utility>
#include <vector>

#include "timur_a_cannon/common/include/common.hpp"
#include "util/include/util.hpp"

namespace timur_a_cannon {

namespace {

using Matrix = std::vector<std::vector<double>>;
using BlockGrid = std::vector<std::vector<Matrix>>;

template <typename Func>
void ParallelFor(int work_size, const Func &func) {
  if (work_size <= 0) {
    return;
  }

  const int num_threads = std::max(1, std::min(ppc::util::GetNumThreads(), work_size));
  std::vector<std::thread> threads;
  threads.reserve(static_cast<std::size_t>(num_threads));

  for (int thread_id = 0; thread_id < num_threads; ++thread_id) {
    threads.emplace_back([&, thread_id]() {
      for (int index = thread_id; index < work_size; index += num_threads) {
        func(index);
      }
    });
  }

  for (auto &thread : threads) {
    thread.join();
  }
}

void DistributeData(const Matrix &src_a, const Matrix &src_b, BlockGrid &bl_a, BlockGrid &bl_b, int b_size,
                    int grid_sz) {
  ParallelFor(grid_sz, [&](int i) {
    for (int j = 0; j < grid_sz; ++j) {
      const int shift = (i + j) % grid_sz;
      for (int row = 0; row < b_size; ++row) {
        for (int col = 0; col < b_size; ++col) {
          bl_a[i][j][row][col] = src_a[(i * b_size) + row][(shift * b_size) + col];
          bl_b[i][j][row][col] = src_b[(shift * b_size) + row][(j * b_size) + col];
        }
      }
    }
  });
}

void RotateBlocksA(BlockGrid &blocks, int grid_sz) {
  ParallelFor(grid_sz, [&](int i) {
    Matrix first_block = std::move(blocks[i][0]);
    for (int j = 0; j < grid_sz - 1; ++j) {
      blocks[i][j] = std::move(blocks[i][j + 1]);
    }
    blocks[i][grid_sz - 1] = std::move(first_block);
  });
}

void RotateBlocksB(BlockGrid &blocks, int grid_sz) {
  ParallelFor(grid_sz, [&](int j) {
    Matrix first_block = std::move(blocks[0][j]);
    for (int i = 0; i < grid_sz - 1; ++i) {
      blocks[i][j] = std::move(blocks[i + 1][j]);
    }
    blocks[grid_sz - 1][j] = std::move(first_block);
  });
}

void CollectResult(const BlockGrid &bl_c, Matrix &result, int b_size, int grid_sz) {
  ParallelFor(grid_sz, [&](int i) {
    for (int j = 0; j < grid_sz; ++j) {
      for (int row = 0; row < b_size; ++row) {
        for (int col = 0; col < b_size; ++col) {
          result[(i * b_size) + row][(j * b_size) + col] = bl_c[i][j][row][col];
        }
      }
    }
  });
}

}  // namespace

TimurACannonMatrixMultiplicationSTL::TimurACannonMatrixMultiplicationSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool TimurACannonMatrixMultiplicationSTL::ValidationImpl() {
  const auto &input = GetInput();
  int b_size = std::get<0>(input);
  const auto &mat_a = std::get<1>(input);
  const auto &mat_b = std::get<2>(input);

  if (b_size <= 0 || mat_a.empty() || mat_b.empty()) {
    return false;
  }

  const std::size_t n = mat_a.size();
  if (mat_b.size() != n || (n % static_cast<std::size_t>(b_size) != 0)) {
    return false;
  }

  const auto is_square_n = [n](const Matrix &matrix) {
    return std::ranges::all_of(matrix, [n](const std::vector<double> &row) { return row.size() == n; });
  };

  return is_square_n(mat_a) && is_square_n(mat_b);
}

bool TimurACannonMatrixMultiplicationSTL::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}

void TimurACannonMatrixMultiplicationSTL::BlockMultiplyAccumulate(const std::vector<std::vector<double>> &a,
                                                                  const std::vector<std::vector<double>> &b,
                                                                  std::vector<std::vector<double>> &c, int b_size) {
  for (int i = 0; i < b_size; ++i) {
    for (int k = 0; k < b_size; ++k) {
      double temp = a[i][k];
      for (int j = 0; j < b_size; ++j) {
        c[i][j] += temp * b[k][j];
      }
    }
  }
}

bool TimurACannonMatrixMultiplicationSTL::RunImpl() {
  const auto &input = GetInput();
  const int b_size = std::get<0>(input);
  const auto &src_a = std::get<1>(input);
  const auto &src_b = std::get<2>(input);
  const int n = static_cast<int>(src_a.size());
  const int grid_sz = n / b_size;

  BlockGrid bl_a(grid_sz, std::vector<Matrix>(grid_sz, Matrix(b_size, std::vector<double>(b_size))));
  BlockGrid bl_b(grid_sz, std::vector<Matrix>(grid_sz, Matrix(b_size, std::vector<double>(b_size))));
  BlockGrid bl_c(grid_sz, std::vector<Matrix>(grid_sz, Matrix(b_size, std::vector<double>(b_size, 0.0))));

  DistributeData(src_a, src_b, bl_a, bl_b, b_size, grid_sz);

  for (int step = 0; step < grid_sz; ++step) {
    ParallelFor(grid_sz, [&](int i) {
      for (int j = 0; j < grid_sz; ++j) {
        BlockMultiplyAccumulate(bl_a[i][j], bl_b[i][j], bl_c[i][j], b_size);
      }
    });

    if (step < grid_sz - 1) {
      RotateBlocksA(bl_a, grid_sz);
      RotateBlocksB(bl_b, grid_sz);
    }
  }

  Matrix result(n, std::vector<double>(n));
  CollectResult(bl_c, result, b_size, grid_sz);

  GetOutput() = std::move(result);
  return true;
}

bool TimurACannonMatrixMultiplicationSTL::PostProcessingImpl() {
  return true;
}

}  // namespace timur_a_cannon
