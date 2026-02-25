#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "buzulukski_d_gaus_gorizontal/common/include/common.hpp"
#include "buzulukski_d_gaus_gorizontal/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"

namespace buzulukski_d_gaus_gorizontal {

class BuzulukskiDGausGorizontalFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<1>(test_param) + "_" + std::to_string(std::get<0>(test_param));
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
  }
  bool CheckTestOutputData(OutType &output_data) final {
    return output_data >= 0 && output_data <= 255;
  }
  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_ = 0;
};

namespace {
bool RunFullPipeline(BuzulukskiDGausGorizontalSEQ &task) {
  return task.Validation() && task.PreProcessing() && task.Run() && task.PostProcessing();
}
}  // namespace

TEST_P(BuzulukskiDGausGorizontalFuncTests, MatmulFromPic) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {std::make_tuple(3, "Small"), std::make_tuple(10, "Medium"),
                                            std::make_tuple(21, "Large")};

const auto kTestTasksList =
    ppc::util::AddFuncTask<BuzulukskiDGausGorizontalSEQ, InType>(kTestParam, PPC_SETTINGS_buzulukski_d_gaus_gorizontal);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);
const auto kTestName = BuzulukskiDGausGorizontalFuncTests::PrintFuncTestName<BuzulukskiDGausGorizontalFuncTests>;

INSTANTIATE_TEST_SUITE_P(Sequential, BuzulukskiDGausGorizontalFuncTests, kGtestValues, kTestName);

TEST(BuzulukskiDGausGorizontalExtra, ConstantImageProcessing) {
  auto task = std::make_shared<BuzulukskiDGausGorizontalSEQ>(10);
  ASSERT_TRUE(RunFullPipeline(*task));
  EXPECT_EQ(task->GetOutput(), 100);
}

TEST(BuzulukskiDGausGorizontalExtra, AllZerosImage) {
  auto task = std::make_shared<BuzulukskiDGausGorizontalSEQ>(5);
  ASSERT_TRUE(task->Validation() && task->PreProcessing());
  std::fill(task->InputImage().begin(), task->InputImage().end(), static_cast<uint8_t>(0));
  ASSERT_TRUE(task->Run() && task->PostProcessing());
  EXPECT_EQ(task->GetOutput(), 0);
}

TEST(BuzulukskiDGausGorizontalExtra, ValidationFailsForSmallSize) {
  auto task = std::make_shared<BuzulukskiDGausGorizontalSEQ>(2);
  EXPECT_FALSE(task->Validation());
}

}  // namespace buzulukski_d_gaus_gorizontal
