#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <tuple>
#include <vector>

#include "makoveeva_matmul_double_seq/common/include/common.hpp"
#include "makoveeva_matmul_double_seq/seq/include/ops_seq.hpp"

namespace makoveeva_matmul_double_seq {
namespace {

// Эталонное умножение матриц - вынесено в отдельную функцию
std::vector<double> ReferenceMultiply(const std::vector<double>& a,
                                      const std::vector<double>& b,
                                      size_t n) {
  std::vector<double> c(n * n, 0.0);
  for (size_t i = 0; i < n; ++i) {
    for (size_t k = 0; k < n; ++k) {
      double tmp = a[(i * n) + k];
      for (size_t j = 0; j < n; ++j) {
        c[(i * n) + j] += tmp * b[(k * n) + j];
      }
    }
  }
  return c;
}

// Проверка результатов с допуском
void CheckResults(const std::vector<double>& actual,
                  const std::vector<double>& expected,
                  size_t n) {
  const double epsilon = 1e-10;
  for (size_t i = 0; i < n * n; ++i) {
    ASSERT_NEAR(actual[i], expected[i], epsilon);
  }
}

}  // namespace

// Тест для матрицы 2x2
TEST(MatmulDoubleFunctionalTest, Multiply2x2) {
  size_t n = 2;

  std::vector<double> a = {1.0, 2.0, 3.0, 4.0};
  std::vector<double> b = {5.0, 6.0, 7.0, 8.0};
  std::vector<double> expected = {19.0, 22.0, 43.0, 50.0};

  auto input = std::make_tuple(n, a, b);
  MatmulDoubleSeqTask task(input);

  ASSERT_TRUE(task.ValidationImpl());
  ASSERT_TRUE(task.PreProcessingImpl());
  ASSERT_TRUE(task.RunImpl());
  ASSERT_TRUE(task.PostProcessingImpl());

  auto result = task.GetResult();
  CheckResults(result, expected, n);
}

// Тест для матрицы 3x3
TEST(MatmulDoubleFunctionalTest, Multiply3x3) {
  size_t n = 3;

  std::vector<double> a = {1.0, 2.0, 3.0,
                           4.0, 5.0, 6.0,
                           7.0, 8.0, 9.0};

  std::vector<double> b = {9.0, 8.0, 7.0,
                           6.0, 5.0, 4.0,
                           3.0, 2.0, 1.0};

  std::vector<double> expected = {30.0, 24.0, 18.0,
                                   84.0, 69.0, 54.0,
                                  138.0, 114.0, 90.0};

  auto input = std::make_tuple(n, a, b);
  MatmulDoubleSeqTask task(input);

  ASSERT_TRUE(task.ValidationImpl());
  ASSERT_TRUE(task.PreProcessingImpl());
  ASSERT_TRUE(task.RunImpl());
  ASSERT_TRUE(task.PostProcessingImpl());

  auto result = task.GetResult();
  CheckResults(result, expected, n);
}

// Тест для матрицы 4x4
TEST(MatmulDoubleFunctionalTest, Multiply4x4) {
  size_t n = 4;

  std::vector<double> a = {1.0,  2.0,  3.0,  4.0,
                           5.0,  6.0,  7.0,  8.0,
                           9.0, 10.0, 11.0, 12.0,
                          13.0, 14.0, 15.0, 16.0};

  std::vector<double> b = {16.0, 15.0, 14.0, 13.0,
                           12.0, 11.0, 10.0,  9.0,
                            8.0,  7.0,  6.0,  5.0,
                            4.0,  3.0,  2.0,  1.0};

  auto expected = ReferenceMultiply(a, b, n);
  auto input = std::make_tuple(n, a, b);
  MatmulDoubleSeqTask task(input);

  ASSERT_TRUE(task.ValidationImpl());
  ASSERT_TRUE(task.PreProcessingImpl());
  ASSERT_TRUE(task.RunImpl());
  ASSERT_TRUE(task.PostProcessingImpl());

  auto result = task.GetResult();
  CheckResults(result, expected, n);
}

}  // namespace makoveeva_matmul_double_seq