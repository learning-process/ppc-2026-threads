#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <tuple>
#include <vector>

#include "sokolov_k_matrix_double_fox_seq/common/include/common.hpp"
#include "sokolov_k_matrix_double_fox_seq/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace sokolov_k_matrix_double_fox_seq {

class SokolovKMatrixDoubleFoxPerfTestsSeq : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kN_ = 400;
  const int kBlock_ = 20;
  InType input_data_;

  void SetUp() override {
    std::vector<double> a(static_cast<std::size_t>(kN_ * kN_), 1.5);
    std::vector<double> b(static_cast<std::size_t>(kN_ * kN_), 2.0);
    input_data_ = std::make_tuple(kN_, kBlock_, a, b);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    const double expected = 1.5 * 2.0 * kN_;
    return std::ranges::all_of(output_data, [expected](double val) { return std::abs(val - expected) <= 1e-6; });
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(SokolovKMatrixDoubleFoxPerfTestsSeq, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, SokolovKMatrixDoubleFoxSEQ>(PPC_SETTINGS_sokolov_k_matrix_double_fox_seq);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = SokolovKMatrixDoubleFoxPerfTestsSeq::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, SokolovKMatrixDoubleFoxPerfTestsSeq, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace sokolov_k_matrix_double_fox_seq
