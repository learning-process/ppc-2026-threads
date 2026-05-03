#include <gtest/gtest.h>

#include <tuple>

#include "nalitov_d_dijkstras_algorithm/common/include/common.hpp"
#include "nalitov_d_dijkstras_algorithm/omp/include/ops_omp.hpp"
#include "nalitov_d_dijkstras_algorithm/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace nalitov_d_dijkstras_algorithm {

class NalitovDDijkstrasAlgorithmPerfTests : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr int kGraphSize = 150;
  InType input_data_{};
  OutType expected_output_{};

  void SetUp() override {
    input_data_ = kGraphSize;
    expected_output_ = input_data_ * (input_data_ - 1) / 2;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return expected_output_ == output_data;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(NalitovDDijkstrasAlgorithmPerfTests, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kSeqPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, NalitovDDijkstrasAlgorithmSeq>(PPC_SETTINGS_nalitov_d_dijkstras_algorithm);
const auto kOmpPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, NalitovDDijkstrasAlgorithmOmp>(PPC_SETTINGS_nalitov_d_dijkstras_algorithm);
const auto kAllPerfTasks = std::tuple_cat(kSeqPerfTasks, kOmpPerfTasks);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = NalitovDDijkstrasAlgorithmPerfTests::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, NalitovDDijkstrasAlgorithmPerfTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace nalitov_d_dijkstras_algorithm
