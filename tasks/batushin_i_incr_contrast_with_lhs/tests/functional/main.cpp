#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>

#include "batushin_i_incr_contrast_with_lhs/common/include/common.hpp"
#include "batushin_i_incr_contrast_with_lhs/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace batushin_i_incr_contrast_with_lhs {

class BatushinIRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int case_num = std::get<0>(params);

    switch (case_num) {
      case 1:
        input_data_ = {50, 100, 150, 200};
        expected_output_ = {0, 85, 170, 255};
        break;
      case 2:
        input_data_ = {10, 200, 50, 180, 30};
        expected_output_ = {0, 255, 54, 228, 27};
        break;
      case 3:
        input_data_ = {120, 130, 125, 128};
        expected_output_ = {0, 255, 128, 204};
        break;
      case 4:
        input_data_ = {70, 160, 90, 200, 110};
        expected_output_ = {0, 177, 40, 255, 79};
        break;
      case 5:
        input_data_ = {0, 64, 128, 192, 255};
        expected_output_ = {0, 64, 128, 192, 255};
        break;
      case 6:
        input_data_ = {100, 100, 100};
        expected_output_ = {128, 128, 128};
        break;
      case 7:
        input_data_ = {0, 255, 0, 255};
        expected_output_ = {0, 255, 0, 255};
        break;
      case 8:
        input_data_ = {0, 0, 0};
        expected_output_ = {128, 128, 128};
        break;
      case 9:
        input_data_ = {255, 255};
        expected_output_ = {128, 128};
        break;
      case 10:
        input_data_ = {255};
        expected_output_ = {128};
        break;
      default:
        break;
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != expected_output_.size()) {
      return false;
    }

    for (size_t i = 0; i < output_data.size(); ++i) {
      if (std::abs(static_cast<int>(output_data[i]) - static_cast<int>(expected_output_[i])) > 1) {
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
  OutType expected_output_;
};

namespace {

TEST_P(BatushinIRunFuncTestsThreads, IncreaseContrastTest) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 10> kTestParams = {std::make_tuple(1, "linear"), std::make_tuple(2, "random"),
                                              std::make_tuple(3, "narrow"), std::make_tuple(4, "mixed"),
                                              std::make_tuple(5, "full"),   std::make_tuple(6, "uniform"),
                                              std::make_tuple(7, "bw"),     std::make_tuple(8, "black"),
                                              std::make_tuple(9, "white"),  std::make_tuple(10, "single")};

const auto kTestTasksList =
    ppc::util::AddFuncTask<BatushinITestTaskSEQ, InType>(kTestParams, PPC_SETTINGS_batushin_i_incr_contrast_with_lhs);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = BatushinIRunFuncTestsThreads::PrintFuncTestName<BatushinIRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(IncreaseContrastTests, BatushinIRunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace batushin_i_incr_contrast_with_lhs
