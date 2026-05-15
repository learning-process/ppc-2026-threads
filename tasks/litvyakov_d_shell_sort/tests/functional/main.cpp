#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdio>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include "litvyakov_d_shell_sort/common/include/common.hpp"
#include "litvyakov_d_shell_sort/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace litvyakov_d_shell_sort {

class LitvyakovDShellSortFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 protected:
  void SetUp() override {
    TestType param = std::get<static_cast<int>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int sz = param;
    input_data_ = std::vector<int>(sz, 0);
    std::random_device rd;
    std::mt19937_64 rng(rd());
    std::uniform_int_distribution<int> dist(std::numeric_limits<int>::min(), std::numeric_limits<int>::max());
    for (auto &i : input_data_) {
      i = dist(rng);
    }
    expected_data_ = input_data_;
    std::ranges::sort(expected_data_);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != expected_data_.size()) {
      return false;
    }

    for (std::size_t i = 0; i < output_data.size(); i++) {
      if (output_data[i] != expected_data_[i]) {
        return false;
      }
    }

    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 public:
  static std::string PrintTestParam(const TestType &param) {
    return "Size_" + std::to_string(param);
  }

 private:
  InType input_data_;
  OutType expected_data_;
};

namespace {

TEST_P(LitvyakovDShellSortFuncTests, shellMergeTest) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 4> kTestParam = {1, 10, 100, 1000};

const auto kTestTasksList =
    ppc::util::AddFuncTask<LitvyakovDShellSortSEQ, InType>(kTestParam, PPC_SETTINGS_litvyakov_d_shell_sort);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kTestName = LitvyakovDShellSortFuncTests::PrintFuncTestName<LitvyakovDShellSortFuncTests>;

INSTANTIATE_TEST_SUITE_P(shellSortFuncTests, LitvyakovDShellSortFuncTests, kGtestValues, kTestName);

}  // namespace

}  // namespace litvyakov_d_shell_sort
