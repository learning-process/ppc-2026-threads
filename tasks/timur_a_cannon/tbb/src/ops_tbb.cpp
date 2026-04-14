#include "timur_a_cannon/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <cstddef>
#include <utility>
#include <vector>

#include "timur_a_cannon/common/include/common.hpp"

namespace timur_a_cannon {

namespace {

using Matrix = std::vector<std::vector<double>>;
using BlockGrid = std::vector<std::vector<Matrix>>;

void MultiplyBlocks(const Matrix &a, const Matrix &b, Matrix &c, int b_size) {
  for (int row = 0; row < b_size; ++row) {
    for (int k = 0; k < b_size; ++k) {
      double temp = a[row][k];
      for (int col = 0; col < b_size; ++col) {
        c[row][col] += temp * b[k][col];
      }
    }
  }
}

void DistributeBlocks(const Matrix &a, const Matrix &b, BlockGrid &bl_a, BlockGrid &bl_b, int b_size, int grid_sz) {
  tbb::parallel_for(tbb::blocked_range2d<int>(0, grid_sz, 0, grid_sz), [&](const tbb::blocked_range2d<int> &r) {
    for (int i = r.rows().begin(); i != r.rows().end(); ++i) {
      for (int j = r.cols().begin(); j != r.cols().end(); ++j) {
        int s = (i + j) % grid_sz;
        for (int row = 0; row < b_size; ++row) {
          for (int col = 0; col < b_size; ++col) {
            bl_a[i][j][row][col] = a[(i * b_size) + row][(s * b_size) + col];
            bl_b[i][j][row][col] = b[(s * b_size) + row][(j * b_size) + col];
          }
        }
      }
    }
  });
}

void RotateAll(BlockGrid &bl_a, BlockGrid &bl_b, int grid_sz) {
  tbb::parallel_for(0, grid_sz, [&](int i) {
    Matrix first = std::move(bl_a[i][0]);
    for (int j = 0; j < grid_sz - 1; ++j) {
      bl_a[i][j] = std::move(bl_a[i][j + 1]);
    }
    bl_a[i][grid_sz - 1] = std::move(first);
  });
  tbb::parallel_for(0, grid_sz, [&](int j) {
    Matrix first = std::move(bl_b[0][j]);
    for (int i = 0; i < grid_sz - 1; ++i) {
      bl_b[i][j] = std::move(bl_b[i + 1][j]);
    }
    bl_b[grid_sz - 1][j] = std::move(first);
  });
}

void AssembleResult(const BlockGrid &bl_c, Matrix &res, int b_size, int grid_sz) {
  tbb::parallel_for(tbb::blocked_range2d<int>(0, grid_sz, 0, grid_sz), [&](const tbb::blocked_range2d<int> &r) {
    for (int i = r.rows().begin(); i != r.rows().end(); ++i) {
      for (int j = r.cols().begin(); j != r.cols().end(); ++j) {
        for (int row = 0; row < b_size; ++row) {
          for (int col = 0; col < b_size; ++col) {
            res[(i * b_size) + row][(j * b_size) + col] = bl_c[i][j][row][col];
          }
        }
      }
    }
  });
}
}  // namespace

TimurACannonMatrixMultiplicationTBB::TimurACannonMatrixMultiplicationTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool TimurACannonMatrixMultiplicationTBB::ValidationImpl() {
  const auto &input = GetInput();
  int b_size = std::get<0>(input);
  const auto &mat_a = std::get<1>(input);
  const auto &mat_b = std::get<2>(input);
  if (b_size <= 0 || mat_a.empty() || mat_b.empty()) {
    return false;
  }
  size_t n = mat_a.size();
  return mat_a[0].size() == n && mat_b.size() == n && (n % static_cast<size_t>(b_size) == 0);
}

bool TimurACannonMatrixMultiplicationTBB::PreProcessingImpl() {
  return true;
}

bool TimurACannonMatrixMultiplicationTBB::RunImpl() {
  const auto &input = GetInput();
  int b_size = std::get<0>(input);
  const auto &matrix_a = std::get<1>(input);
  const auto &matrix_b = std::get<2>(input);
  int n = static_cast<int>(matrix_a.size());
  int grid_sz = n / b_size;

  BlockGrid bl_a(grid_sz, std::vector<Matrix>(grid_sz, Matrix(b_size, std::vector<double>(b_size))));
  BlockGrid bl_b(grid_sz, std::vector<Matrix>(grid_sz, Matrix(b_size, std::vector<double>(b_size))));
  BlockGrid bl_c(grid_sz, std::vector<Matrix>(grid_sz, Matrix(b_size, std::vector<double>(b_size, 0.0))));

  DistributeBlocks(matrix_a, matrix_b, bl_a, bl_b, b_size, grid_sz);

  for (int step = 0; step < grid_sz; ++step) {
    tbb::parallel_for(tbb::blocked_range2d<int>(0, grid_sz, 0, grid_sz), [&](const tbb::blocked_range2d<int> &r) {
      for (int i = r.rows().begin(); i != r.rows().end(); ++i) {
        for (int j = r.cols().begin(); j != r.cols().end(); ++j) {
          MultiplyBlocks(bl_a[i][j], bl_b[i][j], bl_c[i][j], b_size);
        }
      }
    });
    if (grid_sz > 1 && step < grid_sz - 1) {
      RotateAll(bl_a, bl_b, grid_sz);
    }
  }

  Matrix result(n, std::vector<double>(n));
  AssembleResult(bl_c, result, b_size, grid_sz);
  GetOutput() = std::move(result);
  return true;
}

bool TimurACannonMatrixMultiplicationTBB::PostProcessingImpl() {
  return true;
}

}  // namespace timur_a_cannon
