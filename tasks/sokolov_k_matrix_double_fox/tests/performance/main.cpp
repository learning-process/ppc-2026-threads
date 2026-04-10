#include <gtest/gtest.h>

#include "sokolov_k_matrix_double_fox/common/include/common.hpp"
#include "sokolov_k_matrix_double_fox/omp/include/ops_omp.hpp"
#include "sokolov_k_matrix_double_fox/seq/include/ops_seq.hpp"
#include "sokolov_k_matrix_double_fox/stl/include/ops_stl.hpp"
#include "sokolov_k_matrix_double_fox/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace sokolov_k_matrix_double_fox {

class SokolovKMatrixDoubleFoxPerfTestsSeq : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kCount_ = 400;
  InType input_data_{};

  void SetUp() override {
    input_data_ = kCount_;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return input_data_ == output_data;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(SokolovKMatrixDoubleFoxPerfTestsSeq, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, SokolovKMatrixDoubleFoxOMP, SokolovKMatrixDoubleFoxSEQ,
                                                       SokolovKMatrixDoubleFoxSTL, SokolovKMatrixDoubleFoxTBB>(
    PPC_SETTINGS_sokolov_k_matrix_double_fox);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = SokolovKMatrixDoubleFoxPerfTestsSeq::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, SokolovKMatrixDoubleFoxPerfTestsSeq, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace sokolov_k_matrix_double_fox
