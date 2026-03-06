#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "savva_d_monte_carlo/common/include/common.hpp"
#include "savva_d_monte_carlo/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace savva_d_monte_carlo {

class SavvaDRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int test_id = std::get<0>(params);

    std::vector<double> lower_bounds;
    std::vector<double> upper_bounds;
    uint64_t num_points = 0;
    std::function<double(const std::vector<double> &)> f;
    double expected = 0.0;

    switch (test_id) {
      // 2D тесты
      case 0: {  // const_2d_unit - f(x,y) = 1
        lower_bounds = {0.0, 0.0};
        upper_bounds = {1.0, 1.0};
        num_points = 1000000;
        f = [](const std::vector<double> &) { return 1.0; };
        expected = 1.0;  // объем = 1*1 = 1
        break;
      }

      case 1: {  // const_2d_large - f(x,y) = 1 на большой области
        lower_bounds = {-10.0, -5.0};
        upper_bounds = {10.0, 5.0};
        num_points = 1000000;
        f = [](const std::vector<double> &) { return 1.0; };
        expected = (20.0) * (10.0);  // (10-(-10)) * (5-(-5)) = 20 * 10 = 200
        abs_tolerance_ = 0.1;
        break;
      }

      case 2: {  // sum_2d_unit - f(x,y) = x + y
        lower_bounds = {0.0, 0.0};
        upper_bounds = {1.0, 1.0};
        num_points = 1000000;
        f = [](const std::vector<double> &x) { return x[0] + x[1]; };
        // ∫∫(x+y)dxdy = ∫(xdx)∫dy + ∫dx∫(ydy) = (1/2)*1 + 1*(1/2) = 1
        expected = 1.0;
        break;
      }

      case 3: {  // sum_2d_shifted - f(x,y) = x + y на [1,3]x[2,4]
        lower_bounds = {1.0, 2.0};
        upper_bounds = {3.0, 4.0};
        num_points = 1000000;
        f = [](const std::vector<double> &x) { return x[0] + x[1]; };
        // ∫∫(x+y)dxdy = ∫xdx * ∫dy + ∫dx * ∫ydy
        double width_x = 2.0;                              // 3-1
        double width_y = 2.0;                              // 4-2
        double int_x = (3.0 * 3.0 - 1.0 * 1.0) / 2.0;      // ∫xdx от 1 до 3 = (9-1)/2 = 4
        double int_y = (4.0 * 4.0 - 2.0 * 2.0) / 2.0;      // ∫ydy от 2 до 4 = (16-4)/2 = 6
        expected = (int_x * width_y) + (width_x * int_y);  // 4*2 + 2*6 = 8 + 12 = 20
        break;
      }

      case 4: {  // product_2d_unit - f(x,y) = x * y
        lower_bounds = {0.0, 0.0};
        upper_bounds = {1.0, 1.0};
        num_points = 1000000;
        f = [](const std::vector<double> &x) { return x[0] * x[1]; };
        expected = 0.25;  // (1/2)*(1/2) = 0.25
        break;
      }

      case 5: {  // sum_sq_2d_unit - f(x,y) = x^2 + y^2
        lower_bounds = {0.0, 0.0};
        upper_bounds = {1.0, 1.0};
        num_points = 1000000;
        f = [](const std::vector<double> &x) { return (x[0] * x[0]) + (x[1] * x[1]); };
        // ∫x^2dx от 0 до 1 = 1/3, то же для y
        // ∫∫(x^2+y^2)dxdy = (1/3)*1 + 1*(1/3) = 2/3 ≈ 0.6667
        expected = 2.0 / 3.0;
        break;
      }

      case 6: {  // sin_cos_2d - f(x,y) = sin(x)*cos(y)
        lower_bounds = {0.0, 0.0};
        upper_bounds = {std::acos(-1.0), std::acos(-1.0) / 2.0};
        num_points = 500000;
        f = [](const std::vector<double> &x) { return std::sin(x[0]) * std::cos(x[1]); };
        expected = 2.0;
        abs_tolerance_ = 0.1;
        rel_tolerance_ = 0.1;
        break;
      }

      // 3D тесты
      case 7: {  // const_3d_unit - f(x,y,z) = 1
        lower_bounds = {0.0, 0.0, 0.0};
        upper_bounds = {1.0, 1.0, 1.0};
        num_points = 1000000;
        f = [](const std::vector<double> &) { return 1.0; };
        expected = 1.0;  // объем = 1*1*1 = 1
        break;
      }

      case 8: {  // sum_3d_unit - f(x,y,z) = x + y + z
        lower_bounds = {0.0, 0.0, 0.0};
        upper_bounds = {1.0, 1.0, 1.0};
        num_points = 1000000;
        f = [](const std::vector<double> &x) { return x[0] + x[1] + x[2]; };
        // Симметрично, интеграл = 3 * (1/2) = 1.5
        expected = 1.5;
        break;
      }

      case 9: {  // product_3d_unit - f(x,y,z) = x * y * z
        lower_bounds = {0.0, 0.0, 0.0};
        upper_bounds = {1.0, 1.0, 1.0};
        num_points = 1000000;
        f = [](const std::vector<double> &x) { return x[0] * x[1] * x[2]; };
        expected = 0.125;  // (1/2)^3 = 0.125
        break;
      }

      case 10: {  // sum_sq_3d_unit - f(x,y,z) = x^2 + y^2 + z^2
        lower_bounds = {0.0, 0.0, 0.0};
        upper_bounds = {1.0, 1.0, 1.0};
        num_points = 1000000;
        f = [](const std::vector<double> &x) { return (x[0] * x[0]) + (x[1] * x[1]) + (x[2] * x[2]); };
        expected = 1.0;  // 3 * (1/3) = 1
        break;
      }

      // 4D тесты
      case 11: {  // const_4d_unit - f(w,x,y,z) = 1
        lower_bounds = {0.0, 0.0, 0.0, 0.0};
        upper_bounds = {1.0, 1.0, 1.0, 1.0};
        num_points = 1000000;
        f = [](const std::vector<double> &) { return 1.0; };
        expected = 1.0;
        break;
      }

      case 12: {  // sum_4d_unit - f(w,x,y,z) = w + x + y + z
        lower_bounds = {0.0, 0.0, 0.0, 0.0};
        upper_bounds = {1.0, 1.0, 1.0, 1.0};
        num_points = 1000000;
        f = [](const std::vector<double> &x) { return x[0] + x[1] + x[2] + x[3]; };
        expected = 2.0;  // 4 * (1/2) = 2
        break;
      }

      // Граничные случаи
      case 13: {  // small_points - мало точек для константы
        lower_bounds = {0.0, 0.0};
        upper_bounds = {1.0, 1.0};
        num_points = 1000;
        f = [](const std::vector<double> &) { return 1.0; };
        expected = 1.0;
        abs_tolerance_ = 0.5;  // Большой допуск
        break;
      }

      case 14: {                    // zero_volume - нулевой объем по одному измерению
        lower_bounds = {1.0, 1.0};  // x_min = x_max = 1.0
        upper_bounds = {1.0, 2.0};
        num_points = 1000000;
        f = [](const std::vector<double> &) { return 1.0; };
        expected = 0.0;  // объем = 0 * 1 = 0
        abs_tolerance_ = 0.1;
        break;
      }

      case 15: {  // negative_region - отрицательные границы
        lower_bounds = {-2.0, -1.0};
        upper_bounds = {-1.0, 1.0};
        num_points = 1000000;
        f = [](const std::vector<double> &x) { return x[0] + x[1]; };
        // Сложное ожидаемое значение, но можно вычислить
        double width_x = 1.0;                                      // -1 - (-2) = 1
        double width_y = 2.0;                                      // 1 - (-1) = 2
        double int_x = ((-1.0) * (-1.0) - (-2.0) * (-2.0)) / 2.0;  // (1 - 4)/2 = -1.5
        double int_y = ((1.0 * 1.0) - (-1.0) * (-1.0)) / 2.0;      // (1 - 1)/2 = 0
        expected = (int_x * width_y) + (width_x * int_y);          // -1.5*2 + 1*0 = -3
        break;
      }

      case 16: {  // large_range - большая область
        lower_bounds = {-100.0, -100.0};
        upper_bounds = {100.0, 100.0};
        num_points = 1000000;
        f = [](const std::vector<double> &) { return 1.0; };
        expected = 40000.0;  // 200 * 200 = 40000
        abs_tolerance_ = 1.0;
        break;
      }

      default:
        throw std::runtime_error("Unknown test id: " + std::to_string(test_id));
    }

    input_data_ = InputData(lower_bounds, upper_bounds, num_points, std::move(f));
    right_output_ = expected;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    // Для малого числа точек используем больший допуск
    if (input_data_.count_points < 10000ULL) {
      return std::abs(output_data - right_output_) <= 0.5;
    }

    double volume = input_data_.Volume();
    if (std::abs(volume) < 1e-10) {
      return std::abs(output_data) <= 0.1;
    }

    if (std::abs(right_output_) < 1e-4) {
      return std::abs(output_data) <= abs_tolerance_;
    }

    return std::abs(output_data - right_output_) / std::abs(right_output_) <= rel_tolerance_;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType right_output_ = 0.0;
  double abs_tolerance_ = 0.05;
  double rel_tolerance_ = 0.05;
};

namespace {

TEST_P(SavvaDRunFuncTestsThreads, MonteCarloTests) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 17> kTestParam = {
    std::make_tuple(0, "const_2d_unit"),    std::make_tuple(1, "const_2d_large"),  std::make_tuple(2, "sum_2d_unit"),
    std::make_tuple(3, "sum_2d_shifted"),   std::make_tuple(4, "product_2d_unit"), std::make_tuple(5, "sum_sq_2d_unit"),
    std::make_tuple(6, "sin_cos_2d"),       std::make_tuple(7, "const_3d_unit"),   std::make_tuple(8, "sum_3d_unit"),
    std::make_tuple(9, "product_3d_unit"),  std::make_tuple(10, "sum_sq_3d_unit"), std::make_tuple(11, "const_4d_unit"),
    std::make_tuple(12, "sum_4d_unit"),     std::make_tuple(13, "small_points"),   std::make_tuple(14, "zero_volume"),
    std::make_tuple(15, "negative_region"), std::make_tuple(16, "large_range")};

const auto kTestTaskList =
    std::tuple_cat(  // ppc::util::AddFuncTask<SavvaDMonteCarloALL, InType>(kTestParam,
                     // PPC_SETTINGS_savva_d_monte_carlo), ppc::util::AddFuncTask<SavvaDMonteCarloOMP,
                     // InType>(kTestParam, PPC_SETTINGS_savva_d_monte_carlo),
        ppc::util::AddFuncTask<SavvaDMonteCarloSEQ, InType>(kTestParam, PPC_SETTINGS_savva_d_monte_carlo));
// ppc::util::AddFuncTask<SavvaDMonteCarloSTL, InType>(kTestParam, PPC_SETTINGS_savva_d_monte_carlo),
// ppc::util::AddFuncTask<SavvaDMonteCarloTBB, InType>(kTestParam, PPC_SETTINGS_savva_d_monte_carlo)

const auto kGtestValues = ppc::util::ExpandToValues(kTestTaskList);

const auto kPerfTestName = SavvaDRunFuncTestsThreads::PrintFuncTestName<SavvaDRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(MonteCarloTests, SavvaDRunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace savva_d_monte_carlo
