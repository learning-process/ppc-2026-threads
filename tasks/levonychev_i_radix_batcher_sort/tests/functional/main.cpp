#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <string>
#include <tuple>

#include "levonychev_i_radix_batcher_sort/common/include/common.hpp"
#include "levonychev_i_radix_batcher_sort/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace levonychev_i_radix_batcher_sort {

class LevonychevIRadixBatcherSortRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int test_id = std::get<0>(params);
    if (test_id == 1) {
      input_data_ = {5};
    } else if (test_id == 2) {
      input_data_ = {170, 45, 75, 90, 2, 24, 802, 66};
    } else if (test_id == 3) {
      input_data_ = {-170, -45, -75, -90, -2, -24, -802, -66};
    } else if (test_id == 4) {
      input_data_ = {-170, 45, 75, -90, 2, 24, -802, 66};
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    for (size_t i = 1; i < output_data.size(); ++i) {
      if (output_data[i - 1] > output_data[i]) {
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

TEST_P(LevonychevIRadixBatcherSortRunFuncTestsThreads, RadixBatcherSortTests) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 4> kTestParam = {std::make_tuple(1, "one_element"), std::make_tuple(2, "only_posiive"),
                                            std::make_tuple(3, "only_negative"), std::make_tuple(4, "mixed")};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<LevonychevIRadixBatcherSortSEQ, InType>(
    kTestParam, PPC_SETTINGS_levonychev_i_radix_batcher_sort));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName =
    LevonychevIRadixBatcherSortRunFuncTestsThreads::PrintFuncTestName<LevonychevIRadixBatcherSortRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(RadixBatcherSortTests, LevonychevIRadixBatcherSortRunFuncTestsThreads, kGtestValues,
                         kPerfTestName);

}  // namespace

}  // namespace levonychev_i_radix_batcher_sort
