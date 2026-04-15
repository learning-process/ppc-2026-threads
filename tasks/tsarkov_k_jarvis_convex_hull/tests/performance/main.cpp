#include <gtest/gtest.h>

#include "tsarkov_k_jarvis_convex_hull/common/include/common.hpp"
#include "tsarkov_k_jarvis_convex_hull/omp/include/ops_omp.hpp"
#include "tsarkov_k_jarvis_convex_hull/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace tsarkov_k_jarvis_convex_hull {

class TsarkovKRunPerfTestSEQ : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    input_data_.clear();

    for (int x_coord = 1; x_coord < 199; ++x_coord) {
      for (int y_coord = 1; y_coord < 199; ++y_coord) {
        input_data_.push_back(Point{.x = x_coord, .y = y_coord});
      }
    }

    input_data_.push_back(Point{.x = 0, .y = 0});
    input_data_.push_back(Point{.x = 200, .y = 0});
    input_data_.push_back(Point{.x = 200, .y = 200});
    input_data_.push_back(Point{.x = 0, .y = 200});
  }

  bool CheckTestOutputData(OutType &output_data) final {
    const OutType expected_output = {
        Point{.x = 0, .y = 0},
        Point{.x = 200, .y = 0},
        Point{.x = 200, .y = 200},
        Point{.x = 0, .y = 200},
    };

    return output_data == expected_output;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

TEST_P(TsarkovKRunPerfTestSEQ, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, TsarkovKJarvisConvexHullOMP, TsarkovKJarvisConvexHullSEQ>(
        PPC_SETTINGS_tsarkov_k_jarvis_convex_hull);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = TsarkovKRunPerfTestSEQ::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, TsarkovKRunPerfTestSEQ, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace tsarkov_k_jarvis_convex_hull
