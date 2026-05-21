#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "makoveeva_matmul_double_all/all/include/ops_all.hpp"
#include "makoveeva_matmul_double_all/common/include/common.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace makoveeva_matmul_double_all {

// ============================================================================
// ОСНОВНЫЕ ФУНКЦИОНАЛЬНЫЕ ТЕСТЫ
// ============================================================================

class MakoveevaALLRunFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    size_t n = std::get<0>(test_param);
    std::string desc = std::get<1>(test_param);
    return desc + "_" + std::to_string(n) + "x" + std::to_string(n);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());

    size_t n = std::get<0>(params);
    size_t size = n * n;

    const double start = 0.5;
    const double step = 0.75;

    std::vector<double> a(size);
    std::vector<double> b(size);

    for (size_t idx = 0; idx < size; ++idx) {
      a[idx] = start + (static_cast<double>(idx) * step);
    }

    for (size_t idx = 0; idx < size; ++idx) {
      b[idx] = start + (static_cast<double>(size - 1 - idx) * step);
    }

    input_data_ = std::make_tuple(n, a, b);

    std::vector<double> expected(size, 0.0);
    ReferenceMultiply(a, b, expected, n);
    expected_output_ = expected;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (expected_output_.size() != output_data.size()) {
      return false;
    }

    const double epsilon = 1e-10;
    for (size_t idx = 0; idx < expected_output_.size(); ++idx) {
      if (std::abs(expected_output_[idx] - output_data[idx]) > epsilon) {
        return false;
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

  static void ReferenceMultiply(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c,
                                size_t n) {
    for (size_t row = 0; row < n; ++row) {
      for (size_t col = 0; col < n; ++col) {
        double sum = 0.0;
        for (size_t k = 0; k < n; ++k) {
          sum += a[(row * n) + k] * b[(k * n) + col];
        }
        c[(row * n) + col] = sum;
      }
    }
  }

 private:
  InType input_data_;
  OutType expected_output_;
};

namespace {

TEST_P(MakoveevaALLRunFuncTests, MatMulFoxAlg) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 12> kTestParams = {
    std::make_tuple(1, "size_1x1"),    std::make_tuple(2, "size_2x2"),    std::make_tuple(3, "size_3x3"),
    std::make_tuple(4, "size_4x4"),    std::make_tuple(5, "size_5x5"),    std::make_tuple(6, "size_6x6"),
    std::make_tuple(7, "size_7x7"),    std::make_tuple(8, "size_8x8"),    std::make_tuple(9, "size_9x9"),
    std::make_tuple(10, "size_10x10"), std::make_tuple(16, "size_16x16"), std::make_tuple(32, "size_32x32")};

const auto kTestTasksList =
    ppc::util::AddFuncTask<MatmulDoubleAllTask, InType>(kTestParams, PPC_SETTINGS_makoveeva_matmul_double_all);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kTestName = MakoveevaALLRunFuncTests::PrintFuncTestName<MakoveevaALLRunFuncTests>;

INSTANTIATE_TEST_SUITE_P(MatMulFoxAlg_Basic, MakoveevaALLRunFuncTests, kGtestValues, kTestName);

// ============================================================================
// ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ ДЛЯ UNIT ТЕСТОВ
// ============================================================================

void ExpectMatrixNear(const std::vector<double> &actual, const std::vector<double> &expected, double epsilon = 1e-10) {
  ASSERT_EQ(actual.size(), expected.size());
  for (size_t i = 0; i < actual.size(); ++i) {
    EXPECT_NEAR(actual[i], expected[i], epsilon);
  }
}

void ValidateNoNanOrInf(const std::vector<double> &matrix) {
  for (double value : matrix) {
    EXPECT_FALSE(std::isnan(value));
    EXPECT_FALSE(std::isinf(value));
  }
}

void CheckMatrixSymmetric(const std::vector<double> &matrix, size_t n, double epsilon = 1e-10) {
  for (size_t row = 0; row < n; ++row) {
    for (size_t col = 0; col < n; ++col) {
      EXPECT_NEAR(matrix[(row * n) + col], matrix[(col * n) + row], epsilon);
    }
  }
}

std::vector<double> CreateDiagonalMatrix(size_t n, double diag_value) {
  std::vector<double> matrix(n * n, 0.0);
  for (size_t i = 0; i < n; ++i) {
    matrix[(i * n) + i] = diag_value;
  }
  return matrix;
}

void CheckDiagonalMatrixResult(const std::vector<double> &matrix, size_t n, double expected_diag,
                               double epsilon = 1e-10) {
  for (size_t row = 0; row < n; ++row) {
    for (size_t col = 0; col < n; ++col) {
      double expected = (row == col) ? expected_diag : 0.0;
      EXPECT_NEAR(matrix[(row * n) + col], expected, epsilon);
    }
  }
}

}  // namespace

// ============================================================================
// ДОПОЛНИТЕЛЬНЫЕ UNIT ТЕСТЫ ДЛЯ УВЕЛИЧЕНИЯ ПОКРЫТИЯ
// ============================================================================

class MatmulDoubleAllUnitTests : public ::testing::Test {
 protected:
  void SetUp() override {
    // Подготовка общих данных для unit тестов
  }
};

// ТЕСТ 1: Валидация с корректными данными
TEST_F(MatmulDoubleAllUnitTests, ValidationWithValidInput) {
  const size_t n = 4;
  const size_t size = n * n;

  std::vector<double> a(size, 1.0);
  std::vector<double> b(size, 2.0);

  auto input = std::make_tuple(n, a, b);
  MatmulDoubleAllTask task(input);

  EXPECT_TRUE(task.ValidationImpl());
}

// ТЕСТ 2: Валидация с неправильным размером A
TEST_F(MatmulDoubleAllUnitTests, ValidationWithInvalidSizeA) {
  const size_t n = 4;

  std::vector<double> a((n * n) - 1, 1.0);
  std::vector<double> b(n * n, 2.0);

  auto input = std::make_tuple(n, a, b);
  MatmulDoubleAllTask task(input);

  EXPECT_FALSE(task.ValidationImpl());
}

// ТЕСТ 3: Валидация с неправильным размером B
TEST_F(MatmulDoubleAllUnitTests, ValidationWithInvalidSizeB) {
  const size_t n = 4;

  std::vector<double> a(n * n, 1.0);
  std::vector<double> b((n * n) + 5, 2.0);

  auto input = std::make_tuple(n, a, b);
  MatmulDoubleAllTask task(input);

  EXPECT_FALSE(task.ValidationImpl());
}

// ТЕСТ 4: Валидация с нулевым размером
TEST_F(MatmulDoubleAllUnitTests, ValidationWithZeroSize) {
  const size_t n = 0;

  std::vector<double> a;
  std::vector<double> b;

  auto input = std::make_tuple(n, a, b);
  MatmulDoubleAllTask task(input);

  EXPECT_FALSE(task.ValidationImpl());
}

// ТЕСТ 5: Preprocessing инициализирует корректно
TEST_F(MatmulDoubleAllUnitTests, PreprocessingInitializesCorrectly) {
  const size_t n = 3;
  const size_t size = n * n;

  std::vector<double> a(size, 1.5);
  std::vector<double> b(size, 2.5);

  auto input = std::make_tuple(n, a, b);
  MatmulDoubleAllTask task(input);

  EXPECT_TRUE(task.PreProcessingImpl());

  const auto &result = task.GetResult();
  EXPECT_EQ(result.size(), size);
  for (size_t idx = 0; idx < size; ++idx) {
    EXPECT_EQ(result[idx], 0.0);
  }
}

// ТЕСТ 6: Простое матричное умножение 2x2
TEST_F(MatmulDoubleAllUnitTests, SimpleMatrixMultiplication2x2) {
  const size_t n = 2;
  std::vector<double> a = {1.0, 2.0, 3.0, 4.0};
  std::vector<double> b = {5.0, 6.0, 7.0, 8.0};
  std::vector<double> expected = {19.0, 22.0, 43.0, 50.0};

  auto input = std::make_tuple(n, a, b);
  MatmulDoubleAllTask task(input);

  task.PreProcessingImpl();
  task.RunImpl();
  task.PostProcessingImpl();

  ExpectMatrixNear(task.GetResult(), expected);
}

// ТЕСТ 7: Умножение на единичную матрицу
TEST_F(MatmulDoubleAllUnitTests, MultiplicationByIdentity) {
  const size_t n = 3;
  const size_t size = n * n;

  std::vector<double> a(size);
  for (size_t idx = 0; idx < size; ++idx) {
    a[idx] = static_cast<double>(idx + 1);
  }

  std::vector<double> identity(size, 0.0);
  for (size_t row = 0; row < n; ++row) {
    identity[(row * n) + row] = 1.0;
  }

  auto input = std::make_tuple(n, a, identity);
  MatmulDoubleAllTask task(input);

  task.PreProcessingImpl();
  task.RunImpl();
  task.PostProcessingImpl();

  ExpectMatrixNear(task.GetResult(), a);
}

// ТЕСТ 8: Умножение нулевой матрицы
TEST_F(MatmulDoubleAllUnitTests, MultiplicationByZeroMatrix) {
  const size_t n = 3;
  const size_t size = n * n;

  std::vector<double> a(size, 2.5);
  std::vector<double> zero_matrix(size, 0.0);
  std::vector<double> expected(size, 0.0);

  auto input = std::make_tuple(n, a, zero_matrix);
  MatmulDoubleAllTask task(input);

  task.PreProcessingImpl();
  task.RunImpl();
  task.PostProcessingImpl();

  ExpectMatrixNear(task.GetResult(), expected);
}

// ТЕСТ 9: Матрица 4x4
TEST_F(MatmulDoubleAllUnitTests, Matrix4x4) {
  const size_t n = 4;
  const size_t size = n * n;

  std::vector<double> a(size, 1.0);
  std::vector<double> b(size, 1.0);
  std::vector<double> expected(size, static_cast<double>(n));

  auto input = std::make_tuple(n, a, b);
  MatmulDoubleAllTask task(input);

  task.PreProcessingImpl();
  task.RunImpl();
  task.PostProcessingImpl();

  ExpectMatrixNear(task.GetResult(), expected);
}

// ТЕСТ 10: Матрица 8x8 (большой размер)
TEST_F(MatmulDoubleAllUnitTests, Matrix8x8Large) {
  const size_t n = 8;
  const size_t size = n * n;

  std::vector<double> a(size);
  std::vector<double> b(size);

  for (size_t idx = 0; idx < size; ++idx) {
    a[idx] = 0.5 + (static_cast<double>(idx) * 0.1);
    b[idx] = 1.0 - (static_cast<double>(idx) * 0.05);
  }

  auto input = std::make_tuple(n, a, b);
  MatmulDoubleAllTask task(input);

  task.PreProcessingImpl();
  task.RunImpl();
  task.PostProcessingImpl();

  const auto &result = task.GetResult();

  EXPECT_EQ(result.size(), size);
  ValidateNoNanOrInf(result);
}

// ТЕСТ 11: Симметричные матрицы
TEST_F(MatmulDoubleAllUnitTests, SymmetricMatrices) {
  const size_t n = 4;
  const size_t size = n * n;

  std::vector<double> a(size);
  for (size_t row = 0; row < n; ++row) {
    for (size_t col = 0; col < n; ++col) {
      a[(row * n) + col] = 1.0 + static_cast<double>((row == col) ? 2 : 0);
    }
  }

  auto input = std::make_tuple(n, a, a);
  MatmulDoubleAllTask task(input);

  task.PreProcessingImpl();
  task.RunImpl();
  task.PostProcessingImpl();

  CheckMatrixSymmetric(task.GetResult(), n);
}

// ТЕСТ 12: Диагональные матрицы
TEST_F(MatmulDoubleAllUnitTests, DiagonalMatrices) {
  const size_t n = 4;

  auto a = CreateDiagonalMatrix(n, 2.0);
  auto b = CreateDiagonalMatrix(n, 3.0);

  auto input = std::make_tuple(n, a, b);
  MatmulDoubleAllTask task(input);

  task.PreProcessingImpl();
  task.RunImpl();
  task.PostProcessingImpl();

  CheckDiagonalMatrixResult(task.GetResult(), n, 6.0);
}

// ТЕСТ 13: Матрицы с отрицательными значениями
TEST_F(MatmulDoubleAllUnitTests, NegativeValues) {
  const size_t n = 3;
  const size_t size = n * n;

  std::vector<double> a(size);
  std::vector<double> b(size);

  for (size_t idx = 0; idx < size; ++idx) {
    a[idx] = -1.0 - (static_cast<double>(idx) * 0.5);
    b[idx] = 2.0 - (static_cast<double>(idx) * 0.3);
  }

  auto input = std::make_tuple(n, a, b);
  MatmulDoubleAllTask task(input);

  task.PreProcessingImpl();
  task.RunImpl();
  task.PostProcessingImpl();

  ValidateNoNanOrInf(task.GetResult());
}

// ТЕСТ 14: PostProcessing сохраняет результат
TEST_F(MatmulDoubleAllUnitTests, PostProcessingPreservesResult) {
  const size_t n = 2;
  const size_t size = n * n;

  std::vector<double> a = {1.0, 2.0, 3.0, 4.0};
  std::vector<double> b = {5.0, 6.0, 7.0, 8.0};

  auto input = std::make_tuple(n, a, b);
  MatmulDoubleAllTask task(input);

  task.PreProcessingImpl();
  task.RunImpl();

  const auto &before_post = task.GetResult();
  std::vector<double> result_before = before_post;

  task.PostProcessingImpl();

  const auto &after_post = task.GetResult();

  for (size_t idx = 0; idx < size; ++idx) {
    EXPECT_EQ(result_before[idx], after_post[idx]);
  }
}

}  // namespace makoveeva_matmul_double_all
