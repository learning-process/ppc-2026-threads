#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

#include "kopilov_d_vertical_gauss_filter/common/include/common.hpp"
#include "kopilov_d_vertical_gauss_filter/omp/include/ops_omp.hpp"
#include "kopilov_d_vertical_gauss_filter/seq/include/ops_seq.hpp"
#include "kopilov_d_vertical_gauss_filter/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace kopilov_d_vertical_gauss_filter {

class GaussFilterPerfTester : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    constexpr int target_width = 8192;
    constexpr int target_height = 8192;

    source_image_.width = target_width;
    source_image_.height = target_height;
    source_image_.data.resize(target_width * target_height);

    std::mt19937 rand_engine(1337);
    std::uniform_int_distribution<int> color_dist(0, 255);

    std::generate(source_image_.data.begin(), source_image_.data.end(),
                  [&]() { return static_cast<uint8_t>(color_dist(rand_engine)); });
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
                                KopilovDVerticalGaussFilterTBB>(PPC_SETTINGS_kopilov_d_vertical_gauss_filter);

INSTANTIATE_TEST_SUITE_P(GaussFilterPerf, GaussFilterPerfTester, ppc::util::TupleToGTestValues(kPerformanceTasks),
                         GaussFilterPerfTester::CustomPerfTestName);

}  // namespace

}  // namespace kopilov_d_vertical_gauss_filter
