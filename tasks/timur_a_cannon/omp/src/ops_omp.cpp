#include "timur_a_cannon/omp/include/ops_omp.hpp"

#include <omp.h>

#include <cstddef>
#include <utility>
#include <vector>

#include "timur_a_cannon/common/include/common.hpp"

namespace timur_a_cannon {

TimurACannonMatrixMultiplicationOMP::TimurACannonMatrixMultiplicationOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool TimurACannonMatrixMultiplicationOMP::ValidationImpl() {
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

bool TimurACannonMatrixMultiplicationOMP::PreProcessingImpl() {
  return true;
}

bool TimurACannonMatrixMultiplicationOMP::RunImpl() {
  const auto &input = GetInput();
  int b_size = std::get<0>(input);
  const auto &matrix_a = std::get<1>(input);
  const auto &matrix_b = std::get<2>(input);
  int n = static_cast<int>(matrix_a.size());
  int grid_sz = n / b_size;

  using Matrix = std::vector<std::vector<double>>;
  using BlockGrid = std::vector<std::vector<Matrix>>;

  BlockGrid bl_a(grid_sz, std::vector<Matrix>(grid_sz, Matrix(b_size, std::vector<double>(b_size))));
  BlockGrid bl_b(grid_sz, std::vector<Matrix>(grid_sz, Matrix(b_size, std::vector<double>(b_size))));
  BlockGrid bl_c(grid_sz, std::vector<Matrix>(grid_sz, Matrix(b_size, std::vector<double>(b_size, 0.0))));

#pragma omp parallel for collapse(2) default(none) shared(matrix_a, matrix_b, bl_a, bl_b, b_size, grid_sz)
  for (int i = 0; i < grid_sz; ++i) {
    for (int j = 0; j < grid_sz; ++j) {
      int s_a = (i + j) % grid_sz;
      int s_b = (i + j) % grid_sz;
      for (int r = 0; r < b_size; ++r) {
        for (int c = 0; c < b_size; ++c) {
          bl_a[i][j][r][c] = matrix_a[i * b_size + r][s_a * b_size + c];
          bl_b[i][j][r][c] = matrix_b[s_b * b_size + r][j * b_size + c];
        }
      }
    }
  }
  for (int step = 0; step < grid_sz; ++step) {
#pragma omp parallel for collapse(2) default(none) shared(bl_a, bl_b, bl_c, b_size, grid_sz)
    for (int i = 0; i < grid_sz; ++i) {
      for (int j = 0; j < grid_sz; ++j) {
        for (int r = 0; r < b_size; ++r) {
          for (int k = 0; k < b_size; ++k) {
            double temp = bl_a[i][j][r][k];
            for (int c = 0; c < b_size; ++c) {
              bl_c[i][j][r][c] += temp * bl_b[i][j][k][c];
            }
          }
        }
      }
    }

    if (grid_sz > 1) {
#pragma omp parallel default(none) shared(bl_a, bl_b, grid_sz)
      {
#pragma omp for nowait
        for (int i = 0; i < grid_sz; ++i) {
          Matrix tmp = std::move(bl_a[i][0]);
          for (int j = 0; j < grid_sz - 1; ++j) {
            bl_a[i][j] = std::move(bl_a[i][j + 1]);
          }
          bl_a[i][grid_sz - 1] = std::move(tmp);
        }
#pragma omp for
        for (int j = 0; j < grid_sz; ++j) {
          Matrix tmp = std::move(bl_b[0][j]);
          for (int i = 0; i < grid_sz - 1; ++i) {
            bl_b[i][j] = std::move(bl_b[i + 1][j]);
          }
          bl_b[grid_sz - 1][j] = std::move(tmp);
        }
      }
    }
  }

  Matrix result(n, std::vector<double>(n));
#pragma omp parallel for collapse(2) default(none) shared(bl_c, result, b_size, grid_sz)
  for (int i = 0; i < grid_sz; ++i) {
    for (int j = 0; j < grid_sz; ++j) {
      for (int r = 0; r < b_size; ++r) {
        for (int c = 0; c < b_size; ++c) {
          result[i * b_size + r][j * b_size + c] = bl_c[i][j][r][c];
        }
      }
    }
  }

  GetOutput() = std::move(result);
  return true;
}

bool TimurACannonMatrixMultiplicationOMP::PostProcessingImpl() {
  return true;
}

}  // namespace timur_a_cannon
