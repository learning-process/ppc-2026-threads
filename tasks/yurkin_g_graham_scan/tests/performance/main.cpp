#include <gtest/gtest.h>

#include <random>
#include <string>
#include <vector>

#include "util/include/perf_test_util.hpp"
#include "util/include/util.hpp"
#include "yurkin_g_graham_scan/common/include/common.hpp"
#include "yurkin_g_graham_scan/seq/include/ops_seq.hpp"

namespace yurkin_g_graham_scan {

class YurkinGGrahamScanPerfTets : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  const int kCount_ = 2000;
  InType input_data_{};

  void SetUp() override {
    std::mt19937_64 rng(123456789);
    std::uniform_real_distribution<double> dist(-1000.0, 1000.0);
    input_data_.clear();
    input_data_.reserve(kCount_);
    for (int i = 0; i < kCount_; ++i) {
      input_data_.push_back({dist(rng), dist(rng)});
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.empty()) {
      return false;
    }
    if (output_data.size() > input_data_.size()) {
      return false;
    }

    auto cross = [](const Point &a, const Point &b, const Point &c) {
      return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    };
    const size_t m = output_data.size();
    if (m < 3) {
      return true;
    }
    for (size_t i = 0; i < m; ++i) {
      const Point &p0 = output_data[i];
      const Point &p1 = output_data[(i + 1) % m];
      const Point &p2 = output_data[(i + 2) % m];
      if (cross(p0, p1, p2) < -1e-12) {
        return false;
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(YurkinGGrahamScanPerfTets, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, YurkinGGrahamScanSEQ>(PPC_SETTINGS_yurkin_g_graham_scan);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = YurkinGGrahamScanPerfTets::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, YurkinGGrahamScanPerfTets, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace yurkin_g_graham_scan
