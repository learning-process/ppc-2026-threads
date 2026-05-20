#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

#include "kopilov_d_vertical_gauss_filter/common/include/common.hpp"
#include "kopilov_d_vertical_gauss_filter/omp/include/ops_omp.hpp"
#include "kopilov_d_vertical_gauss_filter/seq/include/ops_seq.hpp"
#include "kopilov_d_vertical_gauss_filter/stl/include/ops_stl.hpp"
#include "kopilov_d_vertical_gauss_filter/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace kopilov_d_vertical_gauss_filter {

class GaussFilterPerfTester : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    constexpr int kTargetWidth = 8192;
    constexpr int kTargetHeight = 8192;

    source_image_.width = kTargetWidth;
    source_image_.height = kTargetHeight;
    source_image_.data.resize(static_cast<size_t>(kTargetWidth) * static_cast<size_t>(kTargetHeight));

    std::random_device rd;
    std::mt19937 rand_engine(rd());
    std::uniform_int_distribution<int> color_dist(0, 255);

    std::ranges::generate(source_image_.data, [&]() { return static_cast<uint8_t>(color_dist(rand_engine)); });
  }

  bool CheckTestOutputData(OutType &out) final {
    return out.width == source_image_.width && out.height == source_image_.height &&
           out.data.size() == source_image_.data.size();
  }

  InType GetTestInputData() final {
    return source_image_;
  }

 private:
  InType source_image_;
};

TEST_P(GaussFilterPerfTester, MeasurePerformance) {
  ExecuteTest(GetParam());
}

namespace {

const auto kPerformanceTasks =
    ppc::util::MakeAllPerfTasks<InType, KopilovDVerticalGaussFilterSEQ, KopilovDVerticalGaussFilterOMP,
                                KopilovDVerticalGaussFilterSTL, KopilovDVerticalGaussFilterTBB>(
        PPC_SETTINGS_kopilov_d_vertical_gauss_filter);

INSTANTIATE_TEST_SUITE_P(GaussFilterPerf, GaussFilterPerfTester, ppc::util::TupleToGTestValues(kPerformanceTasks),
                         GaussFilterPerfTester::CustomPerfTestName);

}  // namespace

}  // namespace kopilov_d_vertical_gauss_filter
