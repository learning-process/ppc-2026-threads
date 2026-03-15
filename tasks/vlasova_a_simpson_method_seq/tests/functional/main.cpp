#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <string>
#include <tuple>
#include <vector>

#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"
#include "vlasova_a_simpson_method_seq/common/include/common.hpp"
#include "vlasova_a_simpson_method_seq/seq/include/ops_seq.hpp"

namespace vlasova_a_simpson_method_seq {

class VlasovaASimpsonMethodFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<4>(test_param);
  }

 protected:
  void SetUp() override {
    const auto &params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());

    const auto &a = std::get<0>(params);
    const auto &b = std::get<1>(params);
    const auto &n = std::get<2>(params);
    expected_ = std::get<3>(params);
    const std::string &func_name = std::get<4>(params);

    std::function<double(const std::vector<double> &)> func;

    if (func_name.find("constant") != std::string::npos) {
      func = [](const std::vector<double> &) { return 1.0; };
    } else if (func_name.find("linear") != std::string::npos) {
      func = [](const std::vector<double> &x) { return x[0]; };
    } else if (func_name.find("quadratic") != std::string::npos) {
      func = [](const std::vector<double> &x) { return (x[0] * x[0]) + (x.size() > 1 ? x[1] * x[1] : 0.0); };
    } else if (func_name.find("product") != std::string::npos) {
      func = [](const std::vector<double> &x) {
        double prod = 1.0;
        for (double val : x) {
          prod *= val;
        }
        return prod;
      };
    } else {
      func = [](const std::vector<double> &) { return 0.0; };
    }

    input_data_ = SimpsonTask(func, a, b, n);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::abs(output_data - expected_) < 1e-3;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  double expected_ = 0.0;
};

namespace {

const std::array<TestType, 8> kTestCases = {
    // 1D
    std::make_tuple(std::vector<double>{0.0}, std::vector<double>{1.0}, std::vector<int>{10}, 1.0, "constant_1d"),
    std::make_tuple(std::vector<double>{0.0}, std::vector<double>{1.0}, std::vector<int>{10}, 0.5, "linear_1d"),
    std::make_tuple(std::vector<double>{0.0}, std::vector<double>{1.0}, std::vector<int>{10}, 1.0 / 3.0,
                    "quadratic_1d"),

    // 2D
    std::make_tuple(std::vector<double>{0.0, 0.0}, std::vector<double>{1.0, 1.0}, std::vector<int>{10, 10}, 1.0,
                    "constant_2d"),
    std::make_tuple(std::vector<double>{0.0, 0.0}, std::vector<double>{1.0, 1.0}, std::vector<int>{10, 10}, 0.25,
                    "product_2d"),
    std::make_tuple(std::vector<double>{0.0, 0.0}, std::vector<double>{1.0, 1.0}, std::vector<int>{10, 10}, 2.0 / 3.0,
                    "quadratic_2d"),

    // 3D
    std::make_tuple(std::vector<double>{0.0, 0.0, 0.0}, std::vector<double>{1.0, 1.0, 1.0}, std::vector<int>{8, 8, 8},
                    1.0, "constant_3d"),
    std::make_tuple(std::vector<double>{0.0, 0.0, 0.0}, std::vector<double>{1.0, 1.0, 1.0}, std::vector<int>{8, 8, 8},
                    0.125, "product_3d")};

const auto kTestTasksList =
    ppc::util::AddFuncTask<VlasovaASimpsonMethodSEQ, InType>(kTestCases, PPC_SETTINGS_vlasova_a_simpson_method_seq);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kTestName = VlasovaASimpsonMethodFuncTests::PrintFuncTestName<VlasovaASimpsonMethodFuncTests>;

INSTANTIATE_TEST_SUITE_P(SimpsonMethodTests, VlasovaASimpsonMethodFuncTests, kGtestValues, kTestName);

TEST_P(VlasovaASimpsonMethodFuncTests, GetIntegral) {
  ExecuteTest(GetParam());
}

}  // namespace

}  // namespace vlasova_a_simpson_method_seq
