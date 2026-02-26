#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "bortsova_a_integrals_rectangle_seq/common/include/common.hpp"
#include "bortsova_a_integrals_rectangle_seq/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace bortsova_a_integrals_rectangle_seq {

class BortsovaAIntegralsRectangleFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int test_id = std::get<0>(params);

    switch (test_id) {
      case 0:
        input_data_.func = [](const std::vector<double> &) { return 1.0; };
        input_data_.lower_bounds = {0.0};
        input_data_.upper_bounds = {1.0};
        input_data_.num_steps = 100;
        expected_ = 1.0;
        break;
      case 1:
        input_data_.func = [](const std::vector<double> &x) { return x[0]; };
        input_data_.lower_bounds = {0.0};
        input_data_.upper_bounds = {1.0};
        input_data_.num_steps = 100;
        expected_ = 0.5;
        break;
      case 2:
        input_data_.func = [](const std::vector<double> &x) { return x[0] * x[0]; };
        input_data_.lower_bounds = {0.0};
        input_data_.upper_bounds = {1.0};
        input_data_.num_steps = 1000;
        expected_ = 1.0 / 3.0;
        break;
      case 3:
        input_data_.func = [](const std::vector<double> &x) { return x[0] * x[1]; };
        input_data_.lower_bounds = {0.0, 0.0};
        input_data_.upper_bounds = {1.0, 1.0};
        input_data_.num_steps = 100;
        expected_ = 0.25;
        break;
      case 4:
        input_data_.func = [](const std::vector<double> &x) { return x[0] + x[1]; };
        input_data_.lower_bounds = {0.0, 0.0};
        input_data_.upper_bounds = {1.0, 1.0};
        input_data_.num_steps = 100;
        expected_ = 1.0;
        break;
      case 5:
        input_data_.func = [](const std::vector<double> &) { return 1.0; };
        input_data_.lower_bounds = {0.0, 0.0, 0.0};
        input_data_.upper_bounds = {1.0, 1.0, 1.0};
        input_data_.num_steps = 20;
        expected_ = 1.0;
        break;
      case 6:
        input_data_.func = [](const std::vector<double> &x) { return x[0] + x[1] + x[2]; };
        input_data_.lower_bounds = {0.0, 0.0, 0.0};
        input_data_.upper_bounds = {1.0, 1.0, 1.0};
        input_data_.num_steps = 20;
        expected_ = 1.5;
        break;
      case 7:
        input_data_.func = [](const std::vector<double> &x) { return x[0] * x[0]; };
        input_data_.lower_bounds = {-1.0};
        input_data_.upper_bounds = {1.0};
        input_data_.num_steps = 1000;
        expected_ = 2.0 / 3.0;
        break;
      case 8:
        input_data_.func = [](const std::vector<double> &x) { return x[0] * x[0] + x[1] * x[1]; };
        input_data_.lower_bounds = {0.0, 0.0};
        input_data_.upper_bounds = {1.0, 1.0};
        input_data_.num_steps = 100;
        expected_ = 2.0 / 3.0;
        break;
      case 9:
        input_data_.func = [](const std::vector<double> &) { return 5.0; };
        input_data_.lower_bounds = {0.0};
        input_data_.upper_bounds = {2.0};
        input_data_.num_steps = 1;
        expected_ = 10.0;
        break;
      default:
        break;
    }
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

TEST_P(BortsovaAIntegralsRectangleFuncTests, IntegralComputation) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 10> kTestParam = {
    std::make_tuple(0, "const_1d"),    std::make_tuple(1, "linear_1d"),     std::make_tuple(2, "quadratic_1d"),
    std::make_tuple(3, "product_2d"),  std::make_tuple(4, "sum_2d"),        std::make_tuple(5, "const_3d"),
    std::make_tuple(6, "sum_3d"),      std::make_tuple(7, "neg_bounds_1d"), std::make_tuple(8, "sum_squares_2d"),
    std::make_tuple(9, "single_step"),
};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<BortsovaAIntegralsRectangleSEQ, InType>(
    kTestParam, PPC_SETTINGS_bortsova_a_integrals_rectangle_seq));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kTestName = BortsovaAIntegralsRectangleFuncTests::PrintFuncTestName<BortsovaAIntegralsRectangleFuncTests>;

INSTANTIATE_TEST_SUITE_P(IntegralTests, BortsovaAIntegralsRectangleFuncTests, kGtestValues, kTestName);

TEST(bortsova_a_integrals_rectangle_seq_Validation, NullFunction) {
  IntegralInput input;
  input.lower_bounds = {0.0};
  input.upper_bounds = {1.0};
  input.num_steps = 100;

  BortsovaAIntegralsRectangleSEQ task(input);
  EXPECT_FALSE(task.Validation());

  task.GetInput().func = [](const std::vector<double> &) { return 1.0; };
  task.PreProcessing();
  task.Run();
  task.PostProcessing();
}

TEST(bortsova_a_integrals_rectangle_seq_Validation, EmptyBounds) {
  IntegralInput input;
  input.func = [](const std::vector<double> &) { return 1.0; };
  input.num_steps = 100;

  BortsovaAIntegralsRectangleSEQ task(input);
  EXPECT_FALSE(task.Validation());

  task.GetInput().lower_bounds = {0.0};
  task.GetInput().upper_bounds = {1.0};
  task.PreProcessing();
  task.Run();
  task.PostProcessing();
}

TEST(bortsova_a_integrals_rectangle_seq_Validation, MismatchedBounds) {
  IntegralInput input;
  input.func = [](const std::vector<double> &) { return 1.0; };
  input.lower_bounds = {0.0, 0.0};
  input.upper_bounds = {1.0};
  input.num_steps = 100;

  BortsovaAIntegralsRectangleSEQ task(input);
  EXPECT_FALSE(task.Validation());

  task.GetInput().upper_bounds = {1.0, 1.0};
  task.PreProcessing();
  task.Run();
  task.PostProcessing();
}

TEST(bortsova_a_integrals_rectangle_seq_Validation, ZeroSteps) {
  IntegralInput input;
  input.func = [](const std::vector<double> &) { return 1.0; };
  input.lower_bounds = {0.0};
  input.upper_bounds = {1.0};
  input.num_steps = 0;

  BortsovaAIntegralsRectangleSEQ task(input);
  EXPECT_FALSE(task.Validation());

  task.GetInput().num_steps = 100;
  task.PreProcessing();
  task.Run();
  task.PostProcessing();
}

TEST(bortsova_a_integrals_rectangle_seq_Validation, NegativeSteps) {
  IntegralInput input;
  input.func = [](const std::vector<double> &) { return 1.0; };
  input.lower_bounds = {0.0};
  input.upper_bounds = {1.0};
  input.num_steps = -5;

  BortsovaAIntegralsRectangleSEQ task(input);
  EXPECT_FALSE(task.Validation());

  task.GetInput().num_steps = 100;
  task.PreProcessing();
  task.Run();
  task.PostProcessing();
}

}  // namespace

}  // namespace bortsova_a_integrals_rectangle_seq
