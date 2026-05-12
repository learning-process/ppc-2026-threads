#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <string>
#include <tuple>

#include "kiselev_i_trapezoidal_method_for_multidimensional_integrals/common/include/common.hpp"
#include "kiselev_i_trapezoidal_method_for_multidimensional_integrals/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace kiselev_i_trapezoidal_method_for_multidimensional_integrals {

constexpr double kPi = std::numbers::pi;

class KiselevIRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
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
        input_data_.left_bounds = {0, 0};
        input_data_.right_bounds = {1, 1};
        input_data_.step_n_size = {150, 150};
        input_data_.type_function = 0;
        expected_value = 2.0 / 3.0;
        break;

      case 1:
        input_data_.left_bounds = {0, 0};
        input_data_.right_bounds = {2, 2};
        input_data_.step_n_size = {200, 200};
        input_data_.type_function = 0;
        expected_value = 32.0 / 3.0;
        break;

      case 2:
        input_data_.left_bounds = {-1, -1};
        input_data_.right_bounds = {1, 1};
        input_data_.step_n_size = {200, 200};
        input_data_.type_function = 0;
        expected_value = 8.0 / 3.0;
        break;

      case 3:
        input_data_.left_bounds = {0, 0};
        input_data_.right_bounds = {1, 2};
        input_data_.step_n_size = {200, 200};
        input_data_.type_function = 0;
        expected_value = 10.0 / 3.0;
        break;

      case 4:
        input_data_.left_bounds = {1, 1};
        input_data_.right_bounds = {2, 2};
        input_data_.step_n_size = {150, 150};
        input_data_.type_function = 0;
        expected_value = (14.0 / 3.0);
        break;

      case 5:
        input_data_.left_bounds = {0, 0};
        input_data_.right_bounds = {1, 1};
        input_data_.step_n_size = {150, 150};
        input_data_.type_function = 4;
        expected_value = 1.0;
        break;

      case 6:
        input_data_.left_bounds = {0, 0};
        input_data_.right_bounds = {2, 2};
        input_data_.step_n_size = {150, 150};
        input_data_.type_function = 4;
        expected_value = 8.0;
        break;

      case 7:
        input_data_.left_bounds = {-1, -1};
        input_data_.right_bounds = {1, 1};
        input_data_.step_n_size = {200, 200};
        input_data_.type_function = 4;
        expected_value = 0.0;
        break;

      case 8:
        input_data_.left_bounds = {0, 0};
        input_data_.right_bounds = {1, 2};
        input_data_.step_n_size = {150, 150};
        input_data_.type_function = 4;
        expected_value = 3.0;
        break;

      case 9:
        input_data_.left_bounds = {0, 0};
        input_data_.right_bounds = {kPi / 2, kPi / 2};
        input_data_.step_n_size = {300, 300};
        input_data_.type_function = 2;
        expected_value = (((kPi / 2) * (1 - 0)) + ((kPi / 2) * (1 - 0)));
        break;

      case 10:
        input_data_.left_bounds = {0, 0};
        input_data_.right_bounds = {1, 1};
        input_data_.step_n_size = {200, 200};
        input_data_.type_function = 2;
        expected_value = (((1) * (std::cos(0) - std::cos(1))) + ((1) * (std::sin(1) - std::sin(0))));
        break;

      case 11:
        input_data_.left_bounds = {-1, -1};
        input_data_.right_bounds = {1, 1};
        input_data_.step_n_size = {300, 300};
        input_data_.type_function = 2;
        expected_value = (((2) * (std::cos(-1) - std::cos(1))) + ((2) * (std::sin(1) - std::sin(-1))));
        break;

      case 12:
        input_data_.left_bounds = {0, 0};
        input_data_.right_bounds = {1, 1};
        input_data_.step_n_size = {200, 200};
        input_data_.type_function = 3;
        expected_value = (std::numbers::e - 1) * (std::numbers::e - 1);
        break;

      case 13:
        input_data_.left_bounds = {0, 0};
        input_data_.right_bounds = {2, 2};
        input_data_.step_n_size = {200, 200};
        input_data_.type_function = 3;
        expected_value = (std::exp(2.0) - 1) * (std::exp(2.0) - 1);
        break;

      case 14:
        input_data_.left_bounds = {0, 0};
        input_data_.right_bounds = {0, 0};
        input_data_.step_n_size = {10, 10};
        input_data_.type_function = 3;
        expected_value = 0.0;
        break;

      default:
        throw std::runtime_error("Unknown test id");
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::abs(output_data - expected_value) < 1e-2;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  // double expected_value;
};

TEST_P(KiselevIRunFuncTestsThreads, IntegralCorrectness) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 15> kTestParam = {
    std::make_tuple(0, "sq1"),   std::make_tuple(1, "sq2"),   std::make_tuple(2, "sq3"),   std::make_tuple(3, "sq4"),
    std::make_tuple(4, "sq5"),   std::make_tuple(5, "lin1"),  std::make_tuple(6, "lin2"),  std::make_tuple(7, "lin3"),
    std::make_tuple(8, "lin4"),  std::make_tuple(9, "cos1"),  std::make_tuple(10, "cos2"), std::make_tuple(11, "cos3"),
    std::make_tuple(12, "exp1"), std::make_tuple(13, "exp2"), std::make_tuple(14, "exp3")};

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<KiselevITestTaskSEQ, InType>(kTestParam, PPC_SETTINGS_example_threads));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

namespace {
INSTANTIATE_TEST_SUITE_P(IntegralTests, KiselevIRunFuncTestsThreads, kGtestValues,
                         KiselevIRunFuncTestsThreads::PrintFuncTestName<KiselevIRunFuncTestsThreads>);
}  // namespace
}  // namespace kiselev_i_trapezoidal_method_for_multidimensional_integrals
