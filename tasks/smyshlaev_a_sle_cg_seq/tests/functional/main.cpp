#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>

#include "smyshlaev_a_sle_cg_seq/common/include/common.hpp"
#include "smyshlaev_a_sle_cg_seq/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace smyshlaev_a_sle_cg_seq {

class SmyshlaevASleCgFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int test_id = std::get<0>(params);

    if (test_id == 1) {
      input_data_.A = {{4.0, 1.0}, {1.0, 3.0}};
      input_data_.b = {9.0, 5.0};
      expected_x_ = {2.0, 1.0};
    } else if (test_id == 2) {
      input_data_.A = {{10.0, -1.0, 2.0}, {-1.0, 11.0, -1.0}, {2.0, -1.0, 10.0}};
      input_data_.b = {14.0, 18.0, 30.0};
      expected_x_ = {1.0, 2.0, 3.0};
    } else if (test_id == 3) {
      input_data_.A = {{1.0, 0.0, 0.0, 0.0}, {0.0, 1.0, 0.0, 0.0}, {0.0, 0.0, 1.0, 0.0}, {0.0, 0.0, 0.0, 1.0}};
      input_data_.b = {5.0, -2.0, 8.0, 1.0};
      expected_x_ = {5.0, -2.0, 8.0, 1.0};
    } else if (test_id == 4) {
      input_data_.A = {{4.0, 1.0}, {1.0, 3.0}};
      input_data_.b = {1.0, 2.0};
      expected_x_ = {1.0 / 11.0, 7.0 / 11.0};
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != expected_x_.size()) {
      return false;
    }
    const double epsilon = 1e-6;
    for (size_t i = 0; i < output_data.size(); ++i) {
      if (std::fabs(output_data[i] - expected_x_[i]) > epsilon) {
        return false;
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_x_;
};

namespace {

TEST_P(SmyshlaevASleCgFuncTests, SleCgSeqTests) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 4> kTestParam = {std::make_tuple(1, "Simple2x2"), std::make_tuple(2, "Medium3x3"),
                                            std::make_tuple(3, "IdentityMatrix4x4"), std::make_tuple(4, "Harder2x2")};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<SmyshlaevASleCgTaskSEQ, InType>(kTestParam, PPC_SETTINGS_smyshlaev_a_sle_cg_seq));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = SmyshlaevASleCgFuncTests::PrintFuncTestName<SmyshlaevASleCgFuncTests>;

INSTANTIATE_TEST_SUITE_P(SleCgSeqFunctionalTests, SmyshlaevASleCgFuncTests, kGtestValues, kPerfTestName);

}  // namespace
}  // namespace smyshlaev_a_sle_cg_seq
