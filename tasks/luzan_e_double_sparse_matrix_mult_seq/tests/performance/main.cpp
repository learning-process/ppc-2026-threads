#include <gtest/gtest.h>

#include "luzan_e_double_sparse_matrix_mult_seq/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace luzan_e_double_sparse_matrix_mult_seq {

class LuzanEDoubleSparseMatrixMultSeqPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kCount_ = 200;
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

TEST_P(LuzanEDoubleSparseMatrixMultSeqPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, NesterovATestTaskALL, NesterovATestTaskOMP, LuzanEDoubleSparseMatrixMultSeq,
                                NesterovATestTaskSTL, NesterovATestTaskTBB>(PPC_SETTINGS_example_threads);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = LuzanEDoubleSparseMatrixMultSeqPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, LuzanEDoubleSparseMatrixMultSeqPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace luzan_e_double_sparse_matrix_mult_seq
