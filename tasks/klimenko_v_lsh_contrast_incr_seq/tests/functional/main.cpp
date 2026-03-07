#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <string>

#include "klimenko_v_lsh_contrast_incr_seq/common/include/common.hpp"
#include "klimenko_v_lsh_contrast_incr_seq/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace klimenko_v_lsh_contrast_incr_seq {

class KlimenkoVRunFuncTestsLSH : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(test_param);
  }

 protected:
  void SetUp() override {
    TestType size = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());

    input_data_.resize(size);
    for (int i = 0; i < size; ++i) {
      input_data_[i] = 50 + (i % 101);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != input_data_.size()) {
      return false;
    }

    auto minmax = std::ranges::minmax_element(input_data_);
    int min_val = *minmax.min;
    int max_val = *minmax.max;

    if (max_val == min_val) {
      return output_data == input_data_;
    }

    for (size_t i = 0; i < input_data_.size(); ++i) {
      int expected = ((input_data_[i] - min_val) * 255) / (max_val - min_val);

      if (expected != output_data[i]) {
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
};

namespace {

TEST_P(KlimenkoVRunFuncTestsLSH, ContrastStretching) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {16, 256, 1024};

const auto kTestTasksList = ppc::util::AddFuncTask<KlimenkoVLSHContrastIncrSEQ, InType>(
    kTestParam, PPC_SETTINGS_klimenko_v_lsh_contrast_incr_seq);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kTestName = KlimenkoVRunFuncTestsLSH::PrintFuncTestName<KlimenkoVRunFuncTestsLSH>;

INSTANTIATE_TEST_SUITE_P(ContrastIncrTests, KlimenkoVRunFuncTestsLSH, kGtestValues, kTestName);

}  // namespace

}  // namespace klimenko_v_lsh_contrast_incr_seq
