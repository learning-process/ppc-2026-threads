#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

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
  void SetUp() override {}

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data.size() > 0;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_ = {-170, 45, 75, -90, 2, 24, -802, 66};
};

namespace {

TEST_P(LevonychevIRadixBatcherSortRunFuncTestsThreads, RadixBatcherSortTests) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {std::make_tuple(3, "3"), std::make_tuple(5, "5"), std::make_tuple(7, "7")};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<LevonychevIRadixBatcherSortSEQ, InType>(
    kTestParam, PPC_SETTINGS_levonychev_i_radix_batcher_sort));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName =
    LevonychevIRadixBatcherSortRunFuncTestsThreads::PrintFuncTestName<LevonychevIRadixBatcherSortRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(RadixBatcherSortTests, LevonychevIRadixBatcherSortRunFuncTestsThreads, kGtestValues,
                         kPerfTestName);

}  // namespace

}  // namespace levonychev_i_radix_batcher_sort
