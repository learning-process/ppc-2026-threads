#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <tuple>
#include <vector>

#include "makoveeva_matmul_double_omp/common/include/common.hpp"
#include "makoveeva_matmul_double_omp/omp/include/ops_omp.hpp"
#include "util/include/func_test_util.hpp"

namespace makoveeva_matmul_double_omp {

using TestType = makoveeva_matmul_double_seq::TestType;

class MatmulDoubleOMPFuncTest : public ppc::util::BaseRunFuncTests<makoveeva_matmul_double_seq::InType,
                                                                   makoveeva_matmul_double_seq::OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param) + "_omp";
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    size_t n = std::get<0>(params);
    std::string type = std::get<1>(params);

    size_t size = n * n;
    std::vector<double> a(size);
    std::vector<double> b(size);

    // Генерация тестовых данных
    for (size_t i = 0; i < size; ++i) {
      a[i] = static_cast<double>(i + 1);
      b[i] = static_cast<double>(size - i);
    }

    input_data_ = std::make_tuple(n, a, b);

    // Эталонное умножение
    std::vector<double> expected(size, 0.0);
    for (size_t i = 0; i < n; ++i) {
      for (size_t k = 0; k < n; ++k) {
        double aik = a[i * n + k];
        for (size_t j = 0; j < n; ++j) {
          expected[i * n + j] += aik * b[k * n + j];
        }
      }
    }
    expected_output_ = expected;
  }

  bool CheckTestOutputData(makoveeva_matmul_double_seq::OutType &output_data) final {
    if (expected_output_.size() != output_data.size()) {
      return false;
    }

    const double epsilon = 1e-10;
    for (size_t i = 0; i < expected_output_.size(); ++i) {
      if (std::abs(expected_output_[i] - output_data[i]) > epsilon) {
        return false;
      }
    }
    return true;
  }

  makoveeva_matmul_double_seq::InType GetTestInputData() final {
    return input_data_;
  }

 private:
  makoveeva_matmul_double_seq::InType input_data_;
  makoveeva_matmul_double_seq::OutType expected_output_;
};

namespace {

TEST_P(MatmulDoubleOMPFuncTest, FoxMatmul) {
  ExecuteTest(GetParam());
}

// Массив тестовых параметров - используем правильный тип
const std::array<TestType, 10> kTestParams = {
    std::make_tuple(1, "1"), std::make_tuple(2, "2"),  std::make_tuple(3, "3"), std::make_tuple(4, "4"),
    std::make_tuple(5, "5"), std::make_tuple(6, "6"),  std::make_tuple(7, "7"), std::make_tuple(8, "8"),
    std::make_tuple(9, "9"), std::make_tuple(10, "10")};

// Создаем список задач
const auto kTestTasksList = ppc::util::AddFuncTask<MatmulDoubleOMPTask, makoveeva_matmul_double_seq::InType>(
    kTestParams, PPC_SETTINGS_makoveeva_matmul_double_omp);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kTestName = MatmulDoubleOMPFuncTest::PrintFuncTestName<MatmulDoubleOMPFuncTest>;

INSTANTIATE_TEST_SUITE_P(MatmulOMPTests, MatmulDoubleOMPFuncTest, kGtestValues, kTestName);

}  // namespace

}  // namespace makoveeva_matmul_double_omp
