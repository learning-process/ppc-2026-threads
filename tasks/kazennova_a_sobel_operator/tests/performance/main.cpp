#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>

#include "kazennova_a_sobel_operator/common/include/common.hpp"
#include "kazennova_a_sobel_operator/seq/include/ops_seq.hpp"
// #include "kazennova_a_sobel_operator/omp/include/ops_omp.hpp"
#include "util/include/perf_test_util.hpp"

namespace kazennova_a_sobel_operator {

class KazennovaARunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr size_t kImageSize = 2000;
  InType input_data_;

  void SetUp() override {
    input_data_.resize(kImageSize * kImageSize);
    for (size_t i = 0; i < input_data_.size(); ++i) {
      input_data_[i] = static_cast<uint8_t>(i % 256);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return !output_data.empty() && output_data.size() == input_data_.size();
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(KazennovaARunPerfTestThreads, RunPerfTests) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, SobelSeq>(PPC_SETTINGS_kazennova_a_sobel_operator);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = KazennovaARunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunPerfTests, KazennovaARunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kazennova_a_sobel_operator
