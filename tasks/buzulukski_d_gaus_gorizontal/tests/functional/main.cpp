#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstdint>
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
  static std::string PrintTestParam(
      const ::testing::TestParamInfo<ppc::util::FuncTestParam<InType, OutType, TestType>> &info) {
    const auto &params_tuple = std::get<TestType>(info.param);
    return std::get<std::string>(params_tuple);
  }

 protected:
  void SetUp() override {
    const auto &params_tuple = std::get<TestType>(GetParam());
    input_data_ = std::get<int>(params_tuple);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data >= 0;
  }
  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_ = 0;
};

TEST_P(BuzulukskiDGausGorizontalFuncTests, SequentialRun) {
  ExecuteTest(GetParam());
}

namespace {
const std::array<TestType, 2> kTestParam = {std::make_tuple(3, "size_3"), std::make_tuple(10, "size_10")};

INSTANTIATE_TEST_SUITE_P(
    buzulukski_d_gaus_gorizontal_seq, BuzulukskiDGausGorizontalFuncTests,
    ppc::util::ExpandToValues(std::tuple_cat(ppc::util::AddFuncTask<BuzulukskiDGausGorizontalSEQ, InType>(
        kTestParam, PPC_SETTINGS_buzulukski_d_gaus_gorizontal))),
    BuzulukskiDGausGorizontalFuncTests::PrintTestParam);

TEST(BuzulukskiDGausGorizontalExtra, AllZerosImage) {
  auto task = std::make_shared<BuzulukskiDGausGorizontalSEQ>(4);
  ASSERT_TRUE(task->Validation());
  task->PreProcessing();
  std::fill(task->InputImage().begin(), task->InputImage().end(), static_cast<uint8_t>(0));
  task->Run();
  task->PostProcessing();
  EXPECT_EQ(task->GetOutput(), 0);
}
}  // namespace
}  // namespace buzulukski_d_gaus_gorizontal
