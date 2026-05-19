#include <gtest/gtest.h>
#include <array>
#include <cstddef>
#include <string>
#include <tuple>
#include "solonin_v_radix_sort_batcher/all/include/ops_all.hpp"
#include "solonin_v_radix_sort_batcher/common/include/common.hpp"
#include "solonin_v_radix_sort_batcher/omp/include/ops_omp.hpp"
#include "solonin_v_radix_sort_batcher/seq/include/ops_seq.hpp"
#include "solonin_v_radix_sort_batcher/stl/include/ops_stl.hpp"
#include "solonin_v_radix_sort_batcher/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace solonin_v_radix_sort_batcher {

class RadixBatcherFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &p) {
    return std::to_string(std::get<0>(p)) + "_" + std::get<1>(p);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int id = std::get<0>(params);
    if (id == 1) {
      input_data_ = {42};
    } else if (id == 2) {
      input_data_ = {5, 3, 8, 1, 9, 2, 7, 4, 6};
    } else if (id == 3) {
      input_data_ = {-5, -1, -8, -3, -9, -2, -7};
    } else if (id == 4) {
      input_data_ = {10, -3, 7, -1, 0, 5, -9, 2};
    } else if (id == 5) {
      input_data_ = {INT32_MAX, INT32_MIN, 0, -1, 1};
    } else if (id == 6) {
      input_data_ = {4, 4, 4, 4, 4};
    }
  }

  bool CheckTestOutputData(OutType &output) final {
    for (size_t i = 1; i < output.size(); ++i) {
      if (output[i - 1] > output[i]) return false;
    }
    return true;
  }

  InType GetTestInputData() final { return input_data_; }

 private:
  InType input_data_;
};

namespace {

TEST_P(RadixBatcherFuncTests, RadixBatcherSortTests) { ExecuteTest(GetParam()); }

const std::array<TestType, 6> kParams = {
    std::make_tuple(1, "single"),   std::make_tuple(2, "positive"), std::make_tuple(3, "negative"),
    std::make_tuple(4, "mixed"),    std::make_tuple(5, "extremes"), std::make_tuple(6, "duplicates"),
};

const auto kTasks = std::tuple_cat(
    ppc::util::AddFuncTask<RadixSortBatcherSEQ, InType>(kParams, PPC_SETTINGS_solonin_v_radix_sort_batcher),
    ppc::util::AddFuncTask<RadixSortBatcherOMP, InType>(kParams, PPC_SETTINGS_solonin_v_radix_sort_batcher),
    ppc::util::AddFuncTask<RadixSortBatcherTBB, InType>(kParams, PPC_SETTINGS_solonin_v_radix_sort_batcher),
    ppc::util::AddFuncTask<RadixSortBatcherSTL, InType>(kParams, PPC_SETTINGS_solonin_v_radix_sort_batcher),
    ppc::util::AddFuncTask<RadixSortBatcherALL, InType>(kParams, PPC_SETTINGS_solonin_v_radix_sort_batcher));

const auto kGtestValues = ppc::util::ExpandToValues(kTasks);
const auto kTestName = RadixBatcherFuncTests::PrintFuncTestName<RadixBatcherFuncTests>;

INSTANTIATE_TEST_SUITE_P(RadixBatcherSortTests, RadixBatcherFuncTests, kGtestValues, kTestName);

}  // namespace
}  // namespace solonin_v_radix_sort_batcher
