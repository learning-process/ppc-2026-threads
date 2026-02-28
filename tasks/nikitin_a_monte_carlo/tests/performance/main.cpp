#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <vector>

// #include "nikitin_a_monte_carlo/all/include/ops_all.hpp"
#include "nikitin_a_monte_carlo/common/include/common.hpp"
// #include "nikitin_a_monte_carlo/omp/include/ops_omp.hpp"
#include "nikitin_a_monte_carlo/seq/include/ops_seq.hpp"
// #include "nikitin_a_monte_carlo/stl/include/ops_stl.hpp"
// #include "nikitin_a_monte_carlo/tbb/include/ops_tbb.hpp"
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
  double expected_value_;
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
  double expected_value_;
};

// Тест 3: 5D интегрирование сложной функции (квадратичная) - самый тяжелый
class NikitinAMonteCarloQuadratic5DPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  void SetUp() override {
    // 5D интегрирование f = x1^2 + x2^2 на гиперпараллелепипеде [0,2]^5
    // Точное значение: объем * среднее значение
    // Для каждой из 2 переменных: ∫0^2 x^2 dx = 8/3, по остальным 3 измерениям просто длина 2
    // Итого: (8/3 + 8/3) * 2^3 = (16/3) * 8 = 128/3 ≈ 42.6667
    const std::size_t dim = 5;
    const int num_points = 3000000;  // 3 миллиона точек
    
    std::vector<double> lower(dim, 0.0);
    std::vector<double> upper(dim, 2.0);
    
    input_data_ = std::make_tuple(lower, upper, num_points, FunctionType::kQuadratic);
    expected_value_ = 128.0 / 3.0;  // ≈ 42.6667
  }

  InType GetTestInputData() final {
    return input_data_;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    // Допустимая погрешность 3% для 5D сложной функции
    double relative_error = std::abs(output_data - expected_value_) / expected_value_;
    return relative_error <= 0.03;
  }

 private:
  InType input_data_;
  double expected_value_;
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

INSTANTIATE_TEST_SUITE_P(RunPerfConstant3D, NikitinAMonteCarloConstant3DPerfTests, 
                         kConstant3DGtestValues, kConstant3DPerfTestName);

// Тесты для 4D линейной функции
TEST_P(NikitinAMonteCarloLinear4DPerfTests, RunPerfModesLinear4D) {
  ExecuteTest(GetParam());
}

const auto kLinear4DPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, NikitinAMonteCarloSEQ>(PPC_SETTINGS_nikitin_a_monte_carlo);

const auto kLinear4DGtestValues = ppc::util::TupleToGTestValues(kLinear4DPerfTasks);

const auto kLinear4DPerfTestName = NikitinAMonteCarloLinear4DPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunPerfLinear4D, NikitinAMonteCarloLinear4DPerfTests, 
                         kLinear4DGtestValues, kLinear4DPerfTestName);

}  // namespace

}  // namespace nikitin_a_monte_carlo