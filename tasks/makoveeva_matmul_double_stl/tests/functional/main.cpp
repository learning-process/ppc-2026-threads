#include <gtest/gtest.h>

#include <vector>

#include "makoveeva_matmul_double_stl/stl/include/ops_stl.hpp"

namespace makoveeva_matmul_double_stl {
namespace {

void ReferenceMultiply(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c, size_t n) {
  c.assign(n * n, 0.0);
  for (size_t i = 0; i < n; ++i) {
    for (size_t k = 0; k < n; ++k) {
      const double aik = a[(i * n) + k];
      for (size_t j = 0; j < n; ++j) {
        c[(i * n) + j] += aik * b[(k * n) + j];
      }
    }
  }
}

void CheckTaskExecution(MatmulDoubleSTLTask &task) {
  EXPECT_TRUE(task.ValidationImpl());
  EXPECT_TRUE(task.PreProcessingImpl());
  EXPECT_TRUE(task.RunImpl());
  EXPECT_TRUE(task.PostProcessingImpl());
}

void CheckResult(const std::vector<double> &result, const std::vector<double> &expected) {
  ASSERT_EQ(result.size(), expected.size());
  const double epsilon = 1e-10;
  for (size_t i = 0; i < result.size(); ++i) {
    EXPECT_NEAR(result[i], expected[i], epsilon);
  }
}

}  // namespace

TEST(MatmulDoubleSTLTest, Multiply1x1) {
  const size_t n = 1;
  const std::vector<double> a = {2.0};
  const std::vector<double> b = {3.0};
  std::vector<double> expected(n * n, 0.0);
  ReferenceMultiply(a, b, expected, n);

  auto input = std::make_tuple(n, a, b);
  MatmulDoubleSTLTask task(input);

  CheckTaskExecution(task);
  CheckResult(task.GetResult(), expected);
}

TEST(MatmulDoubleSTLTest, Multiply2x2) {
  const size_t n = 2;
  const std::vector<double> a = {1.0, 2.0, 3.0, 4.0};
  const std::vector<double> b = {5.0, 6.0, 7.0, 8.0};
  std::vector<double> expected(n * n, 0.0);
  ReferenceMultiply(a, b, expected, n);

  auto input = std::make_tuple(n, a, b);
  MatmulDoubleSTLTask task(input);

  CheckTaskExecution(task);
  CheckResult(task.GetResult(), expected);
}

TEST(MatmulDoubleSTLTest, Multiply3x3) {
  const size_t n = 3;
  const std::vector<double> a = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
  const std::vector<double> b = {9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0};
  std::vector<double> expected(n * n, 0.0);
  ReferenceMultiply(a, b, expected, n);

  auto input = std::make_tuple(n, a, b);
  MatmulDoubleSTLTask task(input);

  CheckTaskExecution(task);
  CheckResult(task.GetResult(), expected);
}

TEST(MatmulDoubleSTLTest, Multiply4x4) {
  const size_t n = 4;
  const std::vector<double> a = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0};
  const std::vector<double> b = {16.0, 15.0, 14.0, 13.0, 12.0, 11.0, 10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0};
  std::vector<double> expected(n * n, 0.0);
  ReferenceMultiply(a, b, expected, n);

  auto input = std::make_tuple(n, a, b);
  MatmulDoubleSTLTask task(input);

  CheckTaskExecution(task);
  CheckResult(task.GetResult(), expected);
}

TEST(MatmulDoubleSTLTest, Multiply5x5) {
  const size_t n = 5;
  const size_t size = n * n;
  std::vector<double> a(size);
  std::vector<double> b(size);
  std::vector<double> expected(size, 0.0);

  for (size_t i = 0; i < size; ++i) {
    a[i] = static_cast<double>(i + 1);
    b[i] = static_cast<double>(size - i);
  }
  ReferenceMultiply(a, b, expected, n);

  auto input = std::make_tuple(n, a, b);
  MatmulDoubleSTLTask task(input);

  CheckTaskExecution(task);
  CheckResult(task.GetResult(), expected);
}

TEST(MatmulDoubleSTLTest, Multiply6x6) {
  const size_t n = 6;
  const size_t size = n * n;
  std::vector<double> a(size);
  std::vector<double> b(size);
  std::vector<double> expected(size, 0.0);

  for (size_t i = 0; i < size; ++i) {
    a[i] = static_cast<double>(i + 1);
    b[i] = static_cast<double>(size - i);
  }
  ReferenceMultiply(a, b, expected, n);

  auto input = std::make_tuple(n, a, b);
  MatmulDoubleSTLTask task(input);

  CheckTaskExecution(task);
  CheckResult(task.GetResult(), expected);
}

TEST(MatmulDoubleSTLTest, Multiply7x7) {
  const size_t n = 7;
  const size_t size = n * n;
  std::vector<double> a(size);
  std::vector<double> b(size);
  std::vector<double> expected(size, 0.0);

  for (size_t i = 0; i < size; ++i) {
    a[i] = static_cast<double>(i + 1);
    b[i] = static_cast<double>(size - i);
  }
  ReferenceMultiply(a, b, expected, n);

  auto input = std::make_tuple(n, a, b);
  MatmulDoubleSTLTask task(input);

  CheckTaskExecution(task);
  CheckResult(task.GetResult(), expected);
}

TEST(MatmulDoubleSTLTest, Multiply8x8) {
  const size_t n = 8;
  const size_t size = n * n;
  std::vector<double> a(size);
  std::vector<double> b(size);
  std::vector<double> expected(size, 0.0);

  for (size_t i = 0; i < size; ++i) {
    a[i] = static_cast<double>(i + 1);
    b[i] = static_cast<double>(size - i);
  }
  ReferenceMultiply(a, b, expected, n);

  auto input = std::make_tuple(n, a, b);
  MatmulDoubleSTLTask task(input);

  CheckTaskExecution(task);
  CheckResult(task.GetResult(), expected);
}

TEST(MatmulDoubleSTLTest, Multiply16x16) {
  const size_t n = 16;
  const size_t size = n * n;
  std::vector<double> a(size);
  std::vector<double> b(size);
  std::vector<double> expected(size, 0.0);

  for (size_t i = 0; i < size; ++i) {
    a[i] = static_cast<double>(i + 1);
    b[i] = static_cast<double>(size - i);
  }
  ReferenceMultiply(a, b, expected, n);

  auto input = std::make_tuple(n, a, b);
  MatmulDoubleSTLTask task(input);

  CheckTaskExecution(task);
  CheckResult(task.GetResult(), expected);
}

TEST(MatmulDoubleSTLTest, Multiply32x32) {
  const size_t n = 32;
  const size_t size = n * n;
  std::vector<double> a(size);
  std::vector<double> b(size);
  std::vector<double> expected(size, 0.0);

  for (size_t i = 0; i < size; ++i) {
    a[i] = static_cast<double>(i + 1);
    b[i] = static_cast<double>(size - i);
  }
  ReferenceMultiply(a, b, expected, n);

  auto input = std::make_tuple(n, a, b);
  MatmulDoubleSTLTask task(input);

  CheckTaskExecution(task);
  CheckResult(task.GetResult(), expected);
}

}  // namespace makoveeva_matmul_double_stl
