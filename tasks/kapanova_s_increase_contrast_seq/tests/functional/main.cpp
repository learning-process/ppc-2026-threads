
#include <gtest/gtest.h>

#include <cstdint>
#include <tuple>
#include <vector>

#include "kapanova_s_increase_contrast_seq/common/include/common.hpp"
#include "kapanova_s_increase_contrast_seq/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"

namespace kapanova_s_increase_contrast_seq {

class KapanovaSIncreaseContrastFuncTest : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType& test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int size = std::get<0>(params);
    
    input_data_.resize(size * size);
    for (size_t i = 0; i < input_data_.size(); ++i) {
      input_data_[i] = static_cast<uint8_t>(i % 256);
    }
  }

  bool CheckTestOutputData(OutType& output_data) final {
    if (input_data_.empty() || output_data.empty()) {
      return false;
    }
    if (input_data_.size() != output_data.size()) {
      return false;
    }
    
    auto min_in = *std::min_element(input_data_.begin(), input_data_.end());
    auto max_in = *std::max_element(input_data_.begin(), input_data_.end());
    auto min_out = *std::min_element(output_data.begin(), output_data.end());
    auto max_out = *std::max_element(output_data.begin(), output_data.end());
    
    if (min_in == max_in) {
      return std::all_of(output_data.begin(), output_data.end(), 
                        [min_in](uint8_t v) { return v == min_in; });
    }
    
    return (min_out == 0 && max_out == 255);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

namespace {

TEST_P(KapanovaSIncreaseContrastFuncTest, RunTest) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 5> kTestParam = {
    std::make_tuple(4, "4x4"),
    std::make_tuple(8, "8x8"),
    std::make_tuple(16, "16x16"),
    std::make_tuple(32, "32x32"),
    std::make_tuple(64, "64x64")
};

const auto kTestTasksList = ppc::util::AddFuncTask<KapanovaSIncreaseContrastSEQ, InType>(
    kTestParam, PPC_SETTINGS_kapanova_s_increase_contrast_seq);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = KapanovaSIncreaseContrastFuncTest::PrintFuncTestName<KapanovaSIncreaseContrastFuncTest>;

INSTANTIATE_TEST_SUITE_P(IncreaseContrastTests, KapanovaSIncreaseContrastFuncTest, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kapanova_s_increase_contrast_seq