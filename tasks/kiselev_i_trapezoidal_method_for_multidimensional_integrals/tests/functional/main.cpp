#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "kiselev_i_trapezoidal_method_for_multidimensional_integrals/common/include/common.hpp"
#include "kiselev_i_trapezoidal_method_for_multidimensional_integrals/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace kiselev_i_trapezoidal_method_for_multidimensional_integrals {

constexpr double kPi = 3.1415926535897932384626433832795;
class KiselevIRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    const TestType &params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());

    int test_id = std::get<0>(params);

    switch (test_id) {
      case 0:
        input_data_ = {{0, 0}, {1, 1}, {150, 150}, 0};
        expected_value_ = 2.0 / 3.0;
        break;

      case 1:
        input_data_ = {{0, 0}, {2, 2}, {200, 200}, 0};
        expected_value_ = 32.0 / 3.0;
        break;

      case 2:
        input_data_ = {{-1, -1}, {1, 1}, {200, 200}, 0};
        expected_value_ = 8.0 / 3.0;
        break;

      case 3:
        input_data_ = {{0, 0}, {1, 2}, {200, 200}, 0};
        expected_value_ = 10.0 / 3.0;
        break;

      case 4:
        input_data_ = {{1, 1}, {2, 2}, {150, 150}, 0};
        expected_value_ = (14.0 / 3.0);
        break;

      case 5:
        input_data_ = {{0, 0}, {1, 1}, {150, 150}, 4};
        expected_value_ = 1.0;
        break;

      case 6:
        input_data_ = {{0, 0}, {2, 2}, {150, 150}, 4};
        expected_value_ = 8.0;
        break;

      case 7:
        input_data_ = {{-1, -1}, {1, 1}, {200, 200}, 4};
        expected_value_ = 0.0;
        break;

      case 8:
        input_data_ = {{0, 0}, {1, 2}, {150, 150}, 4};
        expected_value_ = 3.0;
        break;

      case 9:
        input_data_ = {{0, 0}, {2 * kPi, 2 * kPi}, {600, 600}, 1};
        expected_value_ = 0.0;
        break;

      case 10:
        input_data_ = {{0, 0}, {kPi / 2, kPi / 2}, {600, 600}, 1};
        expected_value_ = 1.0;
        break;

      case 11:
        input_data_ = {{0, 0}, {kPi, kPi / 2}, {600, 600}, 1};
        expected_value_ = 2.0;
        break;

      case 12:
        input_data_ = {{0, 0}, {kPi / 2, kPi}, {600, 600}, 1};
        expected_value_ = 0.0;
        break;

      case 13:
        input_data_ = {{0, 0}, {kPi, kPi}, {300, 300}, 2};
        expected_value_ = 2 * kPi;
        break;

      case 14:
        input_data_ = {{0, 0}, {kPi / 2, kPi / 2}, {300, 300}, 2};
        expected_value_ = (kPi / 2) * (1 - 0) + (kPi / 2) * (1 - 0);  // π
        break;

      case 15:
        input_data_ = {{0, 0}, {1, 1}, {200, 200}, 2};
        expected_value_ = (1) * (std::cos(0) - std::cos(1)) + (1) * (std::sin(1) - std::sin(0));
        break;

      case 16:
        input_data_ = {{-1, -1}, {1, 1}, {300, 300}, 2};
        expected_value_ = (2) * (std::cos(-1) - std::cos(1)) + (2) * (std::sin(1) - std::sin(-1));
        break;

      case 17:
        input_data_ = {{0, 0}, {1, 1}, {200, 200}, 3};
        expected_value_ = (std::exp(1) - 1) * (std::exp(1) - 1);
        break;

      case 18:
        input_data_ = {{0, 0}, {2, 2}, {200, 200}, 3};
        expected_value_ = (std::exp(2) - 1) * (std::exp(2) - 1);
        break;

      case 19:
        input_data_ = {{0, 0}, {0, 0}, {10, 10}, 3};
        expected_value_ = 0.0;
        break;

      case 20:
        input_data_ = {{0, 0}, {1, 1}, {150, 150}, 0, -1.0};
        expected_value_ = 2.0 / 3.0;
        break;
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::abs(output_data - expected_value_) < 1e-2;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  double expected_value_;
};

TEST_P(KiselevIRunFuncTestsThreads, IntegralCorrectness) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 21> kTestParam = {
    std::make_tuple(0, "sq1"),    std::make_tuple(1, "sq2"),    std::make_tuple(2, "sq3"),
    std::make_tuple(3, "sq4"),    std::make_tuple(4, "sq5"),    std::make_tuple(5, "lin1"),
    std::make_tuple(6, "lin2"),   std::make_tuple(7, "lin3"),   std::make_tuple(8, "lin4"),
    std::make_tuple(9, "trig1"),  std::make_tuple(10, "trig2"), std::make_tuple(11, "trig3"),
    std::make_tuple(12, "trig4"), std::make_tuple(13, "cos1"),  std::make_tuple(14, "cos2"),
    std::make_tuple(15, "cos3"),  std::make_tuple(16, "cos4"),  std::make_tuple(17, "exp1"),
    std::make_tuple(18, "exp2"),  std::make_tuple(19, "exp3"),  std::make_tuple(19, "null_eps_for_perf")};

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<KiselevITestTaskSEQ, InType>(kTestParam, PPC_SETTINGS_example_threads));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

INSTANTIATE_TEST_SUITE_P(IntegralTests, KiselevIRunFuncTestsThreads, kGtestValues,
                         KiselevIRunFuncTestsThreads::PrintFuncTestName<KiselevIRunFuncTestsThreads>);

}  // namespace kiselev_i_trapezoidal_method_for_multidimensional_integrals
