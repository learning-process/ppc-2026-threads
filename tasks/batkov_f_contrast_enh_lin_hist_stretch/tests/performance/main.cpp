#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <tuple>

#include "batkov_f_contrast_enh_lin_hist_stretch/common/include/common.hpp"
#include "batkov_f_contrast_enh_lin_hist_stretch/omp/include/ops_omp.hpp"
#include "batkov_f_contrast_enh_lin_hist_stretch/seq/include/ops_seq.hpp"
#include "batkov_f_contrast_enh_lin_hist_stretch/stl/include/ops_stl.hpp"
#include "batkov_f_contrast_enh_lin_hist_stretch/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace batkov_f_contrast_enh_lin_hist_stretch {

class BatkovFRunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr size_t kImageSize = 5000;
  InType input_data_;

  void SetUp() override {
    input_data_.resize(kImageSize * kImageSize);
    for (size_t row = 0; row < kImageSize; ++row) {
      for (size_t col = 0; col < kImageSize; ++col) {
        uint8_t value = 100 + ((row + col) % 50);
        input_data_[(row * kImageSize) + col] = value;
      }
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    auto [min_it, max_it] = std::ranges::minmax_element(output_data);
    uint8_t min_out = *min_it;
    uint8_t max_out = *max_it;

    if (min_out == max_out) {
      return true;
    }

    return (min_out <= 1) && (max_out >= 254);
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(BatkovFRunPerfTestThreads, RunPerfTests) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = std::tuple_cat(ppc::util::MakeAllPerfTasks<InType, BatkovFContrastEnhLinHistStretchSEQ>(
                                              PPC_SETTINGS_batkov_f_contrast_enh_lin_hist_stretch),
                                          ppc::util::MakeAllPerfTasks<InType, BatkovFContrastEnhLinHistStretchOMP>(
                                              PPC_SETTINGS_batkov_f_contrast_enh_lin_hist_stretch),
                                          ppc::util::MakeAllPerfTasks<InType, BatkovFContrastEnhLinHistStretchTBB>(
                                              PPC_SETTINGS_batkov_f_contrast_enh_lin_hist_stretch),
                                          ppc::util::MakeAllPerfTasks<InType, BatkovFContrastEnhLinHistStretchSTL>(
                                              PPC_SETTINGS_batkov_f_contrast_enh_lin_hist_stretch));

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = BatkovFRunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunPerfTests, BatkovFRunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace batkov_f_contrast_enh_lin_hist_stretch
