#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <string>
#include <tuple>
#include <vector>

#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"
#include "zyuzin_n_multi_integrals_simpson/common/include/common.hpp"
#include "zyuzin_n_multi_integrals_simpson/omp/include/ops_omp.hpp"
#include "zyuzin_n_multi_integrals_simpson/seq/include/ops_seq.hpp"

namespace zyuzin_n_multi_integrals_simpson {

constexpr double kPi = std::numbers::pi;
constexpr double kTolerance = 1e-3;

class ZyuzinNRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<1>(test_param);
  }
  double expected_value = 0.0;

 protected:
  void SetUp() override {
    const TestType &params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int test_id = std::get<0>(params);

    switch (test_id) {
      case 0:
        // 1D интеграл: интегрирование x на [0, 1] = 0.5
        input_data_.lower_bounds = {0.0};
        input_data_.upper_bounds = {1.0};
        input_data_.n_steps = {100};
        input_data_.func = [](const std::vector<double> &p) { return p[0]; };
        expected_value = 0.5;
        break;

      case 1:
        // 1D интеграл: интегрирование x^2 на [0, 1] = 1/3
        input_data_.lower_bounds = {0.0};
        input_data_.upper_bounds = {1.0};
        input_data_.n_steps = {100};
        input_data_.func = [](const std::vector<double> &p) { return p[0] * p[0]; };
        expected_value = 1.0 / 3.0;
        break;

      case 2:
        // 2D интеграл: интегрирование (x + y) на [0,1]x[0,1] = 1
        input_data_.lower_bounds = {0.0, 0.0};
        input_data_.upper_bounds = {1.0, 1.0};
        input_data_.n_steps = {50, 50};
        input_data_.func = [](const std::vector<double> &p) { return p[0] + p[1]; };
        expected_value = 1.0;
        break;

      case 3:
        // 2D интеграл: интегрирование x*y на [0,1]x[0,1] = 0.25
        input_data_.lower_bounds = {0.0, 0.0};
        input_data_.upper_bounds = {1.0, 1.0};
        input_data_.n_steps = {50, 50};
        input_data_.func = [](const std::vector<double> &p) { return p[0] * p[1]; };
        expected_value = 0.25;
        break;

      case 4:
        // 2D интеграл: интегрирование (x^2 + y^2) на [0,1]x[0,1] = 2/3
        input_data_.lower_bounds = {0.0, 0.0};
        input_data_.upper_bounds = {1.0, 1.0};
        input_data_.n_steps = {50, 50};
        input_data_.func = [](const std::vector<double> &p) { return (p[0] * p[0]) + (p[1] * p[1]); };
        expected_value = 2.0 / 3.0;
        break;

      case 5:
        // 2D интеграл: интегрирование константы 1 на [0,2]x[0,3] = 6 (площадь)
        input_data_.lower_bounds = {0.0, 0.0};
        input_data_.upper_bounds = {2.0, 3.0};
        input_data_.n_steps = {20, 30};
        input_data_.func = [](const std::vector<double> &) { return 1.0; };
        expected_value = 6.0;
        break;

      case 6:
        // 3D интеграл: интегрирование константы 1 на [0,1]^3 = 1 (объем)
        input_data_.lower_bounds = {0.0, 0.0, 0.0};
        input_data_.upper_bounds = {1.0, 1.0, 1.0};
        input_data_.n_steps = {20, 20, 20};
        input_data_.func = [](const std::vector<double> &) { return 1.0; };
        expected_value = 1.0;
        break;

      case 7:
        // 3D интеграл: интегрирование x+y+z на [0,1]^3 = 1.5
        input_data_.lower_bounds = {0.0, 0.0, 0.0};
        input_data_.upper_bounds = {1.0, 1.0, 1.0};
        input_data_.n_steps = {20, 20, 20};
        input_data_.func = [](const std::vector<double> &p) { return p[0] + p[1] + p[2]; };
        expected_value = 1.5;
        break;

      case 8:
        // 1D интеграл: интегрирование sin(x) на [0, pi] = 2
        input_data_.lower_bounds = {0.0};
        input_data_.upper_bounds = {kPi};
        input_data_.n_steps = {100};
        input_data_.func = [](const std::vector<double> &p) { return std::sin(p[0]); };
        expected_value = 2.0;
        break;

      case 9:
        // 2D интеграл: интегрирование sin(x)*cos(y) на [0, pi/2]x[0, pi/2] = 1
        input_data_.lower_bounds = {0.0, 0.0};
        input_data_.upper_bounds = {kPi / 2, kPi / 2};
        input_data_.n_steps = {50, 50};
        input_data_.func = [](const std::vector<double> &p) { return std::sin(p[0]) * std::cos(p[1]); };
        expected_value = 1.0;
        break;

      case 10:
        // 2D интеграл: интегрирование exp(x+y) на [0,1]x[0,1] = (e-1)^2
        input_data_.lower_bounds = {0.0, 0.0};
        input_data_.upper_bounds = {1.0, 1.0};
        input_data_.n_steps = {50, 50};
        input_data_.func = [](const std::vector<double> &p) { return std::exp(p[0] + p[1]); };
        expected_value = (std::numbers::e - 1.0) * (std::numbers::e - 1.0);
        break;

      case 11:
        // 3D интеграл: интегрирование x*y*z на [0,1]^3 = 1/8
        input_data_.lower_bounds = {0.0, 0.0, 0.0};
        input_data_.upper_bounds = {1.0, 1.0, 1.0};
        input_data_.n_steps = {20, 20, 20};
        input_data_.func = [](const std::vector<double> &p) { return p[0] * p[1] * p[2]; };
        expected_value = 1.0 / 8.0;
        break;

      default:
        input_data_.lower_bounds = {0.0};
        input_data_.upper_bounds = {1.0};
        input_data_.n_steps = {10};
        input_data_.func = [](const std::vector<double> &) { return 1.0; };
        expected_value = 1.0;
        break;
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::abs(output_data - expected_value) < kTolerance;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

namespace {

TEST_P(ZyuzinNRunFuncTestsThreads, SimpsonMultiDimTest) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 12> kTestParam = {
    std::make_tuple(0, "1d_linear"),   std::make_tuple(1, "1d_quadratic"),   std::make_tuple(2, "2d_sum"),
    std::make_tuple(3, "2d_product"),  std::make_tuple(4, "2d_sum_squares"), std::make_tuple(5, "2d_constant"),
    std::make_tuple(6, "3d_constant"), std::make_tuple(7, "3d_sum"),         std::make_tuple(8, "1d_sin"),
    std::make_tuple(9, "2d_sin_cos"),  std::make_tuple(10, "2d_exp"),        std::make_tuple(11, "3d_product")};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<ZyuzinNSimpsonSEQ, InType>(kTestParam, PPC_SETTINGS_zyuzin_n_multi_integrals_simpson),
    ppc::util::AddFuncTask<ZyuzinNSimpsonOMP, InType>(kTestParam, PPC_SETTINGS_zyuzin_n_multi_integrals_simpson));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = ZyuzinNRunFuncTestsThreads::PrintFuncTestName<ZyuzinNRunFuncTestsThreads>;
INSTANTIATE_TEST_SUITE_P(SimpsonMultiDimTests, ZyuzinNRunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace zyuzin_n_multi_integrals_simpson
