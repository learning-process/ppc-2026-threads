#include <gtest/gtest.h>

#include <random>
#include <vector>

#include "ilin_a_algorithm_graham/common/include/common.hpp"
#include "ilin_a_algorithm_graham/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace ilin_a_algorithm_graham {

class IlinAGrahamPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
 public:
  void SetUp() override {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-100.0, 100.0);

    std::vector<Point> points;
    for (int i = 0; i < 1000; ++i) {
      points.push_back({dis(gen), dis(gen)});
    }

    input_data_ = InputData{.points = std::move(points)};
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return !output_data.empty();
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

TEST_P(IlinAGrahamPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, IlinAGrahamSEQ>(PPC_SETTINGS_ilin_a_algorithm_graham);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = IlinAGrahamPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, IlinAGrahamPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace ilin_a_algorithm_graham
