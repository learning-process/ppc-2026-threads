#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <tuple>
#include <vector>

#include "nikitin_a_monte_carlo/common/include/common.hpp"
#include "nikitin_a_monte_carlo/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace nikitin_a_monte_carlo {

// Тест 1: 3D интегрирование константной функции (самый легкий)
class NikitinAMonteCarloConstant3DPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  void SetUp() override {
    // 3D интегрирование константы f(x,y,z)=1 на кубе [0,10]^3
    // Точное значение интеграла: объем = 10*10*10 = 1000
    const std::size_t dim = 3;
    const int num_points = 1000000;  // 1 миллион точек

    std::vector<double> lower(dim, 0.0);
    std::vector<double> upper(dim, 10.0);

    input_data_ = std::make_tuple(lower, upper, num_points, FunctionType::kConstant);
    expected_value_ = 1000.0;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    // Для perf-тестов проверяем, что результат близок к ожидаемому
    // Допустимая погрешность 1% для 1 млн точек
    double relative_error = std::abs(output_data - expected_value_) / expected_value_;
    return relative_error <= 0.01;
  }

 private:
  InType input_data_;
  double expected_value_ = 0.0;
};

// Тест 2: 4D интегрирование линейной функции (средний по сложности)
class NikitinAMonteCarloLinear4DPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  void SetUp() override {
    // 4D интегрирование f=x1 на гиперпараллелепипеде [0,5]^4
    // Точное значение: объем * среднее значение = (5^4) * (2.5) = 625 * 2.5 = 1562.5
    const std::size_t dim = 4;
    const int num_points = 2000000;  // 2 миллиона точек

    std::vector<double> lower(dim, 0.0);
    std::vector<double> upper(dim, 5.0);

    input_data_ = std::make_tuple(lower, upper, num_points, FunctionType::kLinear);
    expected_value_ = 1562.5;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    // Допустимая погрешность 2% для 4D
    double relative_error = std::abs(output_data - expected_value_) / expected_value_;
    return relative_error <= 0.02;
  }

 private:
  InType input_data_;
  double expected_value_ = 0.0;
};

namespace {

// Тесты для 3D константной функции
TEST_P(NikitinAMonteCarloConstant3DPerfTests, RunPerfModesConstant3D) {
  ExecuteTest(GetParam());
}

const auto kConstant3DPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, NikitinAMonteCarloSEQ>(PPC_SETTINGS_nikitin_a_monte_carlo);

const auto kConstant3DGtestValues = ppc::util::TupleToGTestValues(kConstant3DPerfTasks);

const auto kConstant3DPerfTestName = NikitinAMonteCarloConstant3DPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunPerfConstant3D, NikitinAMonteCarloConstant3DPerfTests, kConstant3DGtestValues,
                         kConstant3DPerfTestName);

// Тесты для 4D линейной функции
TEST_P(NikitinAMonteCarloLinear4DPerfTests, RunPerfModesLinear4D) {
  ExecuteTest(GetParam());
}

const auto kLinear4DPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, NikitinAMonteCarloSEQ>(PPC_SETTINGS_nikitin_a_monte_carlo);

const auto kLinear4DGtestValues = ppc::util::TupleToGTestValues(kLinear4DPerfTasks);

const auto kLinear4DPerfTestName = NikitinAMonteCarloLinear4DPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunPerfLinear4D, NikitinAMonteCarloLinear4DPerfTests, kLinear4DGtestValues,
                         kLinear4DPerfTestName);

}  // namespace

}  // namespace nikitin_a_monte_carlo
