#include <gtest/gtest.h>

#include <cstddef>
#include <random>

#include "klimenko_v_lsh_contrast_incr_seq/common/include/common.hpp"
#include "klimenko_v_lsh_contrast_incr_seq/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace klimenko_v_lsh_contrast_incr_seq {

class KlimenkoVLSHContrastIncrSEQPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  void SetUp() override {
    const size_t size = 10000000;
    input_data_.resize(size);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0, 255);
    for (auto &pixel : input_data_) {
      pixel = dist(gen);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data.size() == input_data_.size();
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

TEST_P(KlimenkoVLSHContrastIncrSEQPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, KlimenkoVLSHContrastIncrSEQ>(PPC_SETTINGS_klimenko_v_lsh_contrast_incr_seq);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = KlimenkoVLSHContrastIncrSEQPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, KlimenkoVLSHContrastIncrSEQPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace klimenko_v_lsh_contrast_incr_seq
