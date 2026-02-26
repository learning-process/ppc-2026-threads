#include <gtest/gtest.h>

#include <cstddef>
#include <cstdlib>
#include <random>
#include <vector>

#include "perepelkin_i_convex_hull_graham_scan/common/include/common.hpp"
#include "perepelkin_i_convex_hull_graham_scan/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace perepelkin_i_convex_hull_graham_scan {

class PerepelkinIConvexHullGrahamScanPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;
  OutType expected_output_;

  size_t base_length_ = std::pow(10, 6);
  size_t scale_factor_ = 2;
  unsigned int seed_ = 42;

  void SetUp() override {
    std::tie(input_data_, expected_output_) = GenerateTestData(base_length_, scale_factor_, seed_);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return expected_output_ == output_data;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

  static std::tuple<InType, OutType> GenerateTestData(size_t base_length, size_t scale_factor, unsigned int seed) {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(-1000.0, 1000.0);

    // Generate a base set of random points
    InType base(base_length);
    for (auto &value : base) {
      value = std::make_pair(dist(gen), dist(gen));
    }

    // Scale up the base set
    InType data;
    data.reserve(base_length * scale_factor);
    for (size_t i = 0; i < scale_factor; ++i) {
      data.insert(data.end(), base.begin(), base.end());
    }

    // Compute the expected output using the sequential implementation
    ppc::task::TaskPtr<InType, OutType> task = std::make_shared<PerepelkinIConvexHullGrahamScanSEQ>(data);
    task->Validation();
    task->PreProcessing();
    task->Run();
    task->PostProcessing();
    OutType expected_output = task->GetOutput();
  
    return {data, expected_output};
  }
};

TEST_P(PerepelkinIConvexHullGrahamScanPerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, PerepelkinIConvexHullGrahamScanSEQ>(PPC_SETTINGS_perepelkin_i_convex_hull_graham_scan);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = PerepelkinIConvexHullGrahamScanPerfTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, PerepelkinIConvexHullGrahamScanPerfTest, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace perepelkin_i_convex_hull_graham_scan
