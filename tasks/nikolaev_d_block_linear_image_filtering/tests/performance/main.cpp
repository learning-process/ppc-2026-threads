#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <tuple>
#include <vector>

#include "nikolaev_d_block_linear_image_filtering/all/include/ops_all.hpp"
#include "nikolaev_d_block_linear_image_filtering/common/include/common.hpp"
#include "nikolaev_d_block_linear_image_filtering/omp/include/ops_omp.hpp"
#include "nikolaev_d_block_linear_image_filtering/seq/include/ops_seq.hpp"
#include "nikolaev_d_block_linear_image_filtering/stl/include/ops_stl.hpp"
#include "nikolaev_d_block_linear_image_filtering/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace nikolaev_d_block_linear_image_filtering {

class NikolaevDBlockLinearImageFilteringPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kWidth_ = 7680;
  const int kHeight_ = 4320;
  InType input_data_;

  void SetUp() override {
    std::vector<uint8_t> in_data(static_cast<size_t>(kWidth_ * kHeight_ * 3), 128);
    input_data_ = std::make_tuple(kWidth_, kHeight_, in_data);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return std::get<2>(input_data_).size() == output_data.size();
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(NikolaevDBlockLinearImageFilteringPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, NikolaevDBlockLinearImageFilteringSEQ, NikolaevDBlockLinearImageFilteringOMP,
                                NikolaevDBlockLinearImageFilteringTBB, NikolaevDBlockLinearImageFilteringSTL,
                                NikolaevDBlockLinearImageFilteringALL>(
        PPC_SETTINGS_nikolaev_d_block_linear_image_filtering);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = NikolaevDBlockLinearImageFilteringPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, NikolaevDBlockLinearImageFilteringPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace nikolaev_d_block_linear_image_filtering
