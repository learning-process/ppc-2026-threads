#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <tuple>
#include <vector>

#include "makoveeva_matmul_double_seq/seq/include/ops_seq.hpp"

namespace makoveeva_matmul_double_seq {
namespace {

void ReferenceMultiply(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c, size_t n) {
  for (size_t i = 0; i < n; ++i) {
    for (size_t k = 0; k < n; ++k) {
      const double tmp = a[(i * n) + k];
      for (size_t j = 0; j < n; ++j) {
        c[(i * n) + j] += tmp * b[(k * n) + j];
      }
    }
  }
}

// Разбиваем функцию на еще более мелкие части
void ValidateTask(MatmulDoubleSeqTask &task) {
  ASSERT_TRUE(task.ValidationImpl());
}

void PreProcessTask(MatmulDoubleSeqTask &task) {
  ASSERT_TRUE(task.PreProcessingImpl());
}

void RunTask(MatmulDoubleSeqTask &task) {
  ASSERT_TRUE(task.RunImpl());
}

void PostProcessTask(MatmulDoubleSeqTask &task) {
  ASSERT_TRUE(task.PostProcessingImpl());
}

void CheckTaskExecution(MatmulDoubleSeqTask &task) {
  ValidateTask(task);
  PreProcessTask(task);
  RunTask(task);
  PostProcessTask(task);
}

void CheckResults(const std::vector<double> &result, const std::vector<double> &expected) {
  ASSERT_EQ(result.size(), expected.size());
  const double epsilon = 1e-10;
  for (size_t i = 0; i < result.size(); ++i) {
    ASSERT_NEAR(result[i], expected[i], epsilon);
  }
}

}  // namespace

TEST(MatmulDoubleFunctionalTest, multiply2x2) {
  const size_t n = 2;
  const std::vector<double> a = {1.0, 2.0, 3.0, 4.0};
  const std::vector<double> b = {5.0, 6.0, 7.0, 8.0};
  const std::vector<double> expected = {19.0, 22.0, 43.0, 50.0};

  auto input = std::make_tuple(n, a, b);
  MatmulDoubleSeqTask task(input);

  CheckTaskExecution(task);
  CheckResults(task.GetResult(), expected);
}

TEST(MatmulDoubleFunctionalTest, multiply3x3) {
  const size_t n = 3;
  const std::vector<double> a = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
  const std::vector<double> b = {9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0};
  const std::vector<double> expected = {30.0, 24.0, 18.0, 84.0, 69.0, 54.0, 138.0, 114.0, 90.0};

  auto input = std::make_tuple(n, a, b);
  MatmulDoubleSeqTask task(input);

  CheckTaskExecution(task);
  CheckResults(task.GetResult(), expected);
}

TEST(MatmulDoubleFunctionalTest, multiply4x4) {
  const size_t n = 4;
  const std::vector<double> a = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0};
  const std::vector<double> b = {16.0, 15.0, 14.0, 13.0, 12.0, 11.0, 10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0};

  std::vector<double> expected(n * n, 0.0);
  ReferenceMultiply(a, b, expected, n);

  auto input = std::make_tuple(n, a, b);
  MatmulDoubleSeqTask task(input);

  CheckTaskExecution(task);
  CheckResults(task.GetResult(), expected);
}

}  // namespace makoveeva_matmul_double_seq
