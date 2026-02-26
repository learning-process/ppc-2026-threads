#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include "shemetov_d_odd_even_mergesort_seq/common/include/common.hpp"
#include "shemetov_d_odd_even_mergesort_seq/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace shemetov_d_odd_even_mergesort {

class ShemetovDRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &param) {
    return "SIZE_" + std::to_string(std::get<0>(param));
  }

 protected:
  void SetUp() override {
    test_array_ = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
  }

  InType GetTestInputData() final {
    return test_array_;
  }

  bool CheckTestOutputData(OutType &output) final {
    return std::ranges::is_sorted(output.begin(), output.end());
  }

 private:
  InType test_array_;
};

namespace {

void ValidateInput(const size_t size, ShemetovDOddEvenMergeSortSEQ &test_task) {
  if (size == 0) {
    ASSERT_FALSE(test_task.Validation());
    return;
  }

  ASSERT_TRUE(test_task.Validation());
}

void PrepareInput(ShemetovDOddEvenMergeSortSEQ &test_task) {
  ASSERT_TRUE(test_task.PreProcessing());
}

void RunTestsuit(ShemetovDOddEvenMergeSortSEQ &test_task) {
  ASSERT_TRUE(test_task.Run());
  ASSERT_TRUE(test_task.PostProcessing());
}

void CheckOutput(const std::vector<int> &array, ShemetovDOddEvenMergeSortSEQ &test_task) {
  OutType output = test_task.GetOutput();

  ASSERT_EQ(output.size(), array.size());
  ASSERT_TRUE(std::ranges::is_sorted(output.begin(), output.end()));
}

void RunOddEvenMergeSortTest(const InType &input) {
  const auto &[size, array] = input;

  ShemetovDOddEvenMergeSortSEQ test_task(input);

  ValidateInput(size, test_task);
  PrepareInput(test_task);
  RunTestsuit(test_task);
  CheckOutput(array, test_task);
}

TEST_P(ShemetovDRunFuncTestsThreads, OddEvenMergeSortFuncTest) {
  RunOddEvenMergeSortTest(std::get<static_cast<size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam()));
}

const std::array<TestType, 12> kTestParam = {
    TestType{0, {}},
    TestType{1, {0}},
    TestType{2, {1, 0}},
    TestType{3, {1, -1, 0}},
    TestType{4, {4, 3, 2, 1}},
    TestType{5, {-12, -12, -12, -12, -12}},
    TestType{7, {0, 6, 1, 5, 2, 4, 3}},
    TestType{8, {10, -1, 3, 0, 0, 4, 5, 9}},
    TestType{11, {-5, 4, 5, -7, 7, 89, 12, 15, 0, 10, 11}},
    TestType{16, {0, 15, 1, 14, 2, 13, 3, 12, 4, 11, 5, 10, 6, 9, 7, 8}},
    TestType{17, {-1, -2, -3, -4, -5, -6, -7, -8, 0, 8, 7, 6, 5, 4, 3, 2, 1}},
    TestType{18, {17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0}}};

const auto kTestTasksList = ppc::util::AddFuncTask<ShemetovDOddEvenMergeSortSEQ, InType>(
    kTestParam, PPC_SETTINGS_shemetov_d_odd_even_mergesort_seq);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = ShemetovDRunFuncTestsThreads::PrintFuncTestName<ShemetovDRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(OddEvenMergeSortTests, ShemetovDRunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace shemetov_d_odd_even_mergesort
