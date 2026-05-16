#include <gtest/gtest.h>

#include <tuple>
#include <vector>

#include "luzan_e_double_sparse_matrix_mult_seq/common/include/common.hpp"
#include "luzan_e_double_sparse_matrix_mult_seq/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace luzan_e_double_sparse_matrix_mult_seq {

class LuzanEDoubleSparseMatrixMultSeqPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kCount_ = 2500;
  // const unsigned cols_ = 5000;
  // const unsigned rows_ = 5000;
  InType input_data_;
  OutType output_data_;

  void SetUp() override {
    std::vector<double> value;
    std::vector<unsigned> row;
    std::vector<unsigned> col_index;

    SparseMatrix lhs;
    lhs.GenLineMatrix(kCount_, kCount_);

    SparseMatrix rhs;
    rhs.GenColsMatrix(kCount_, kCount_);
    input_data_ = std::make_tuple(lhs, rhs);
    output_data_.GenPerfAns(kCount_, kCount_, kCount_);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return (output_data == output_data_);
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(LuzanEDoubleSparseMatrixMultSeqPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, LuzanEDoubleSparseMatrixMultSeq>(
    PPC_SETTINGS_luzan_e_double_sparse_matrix_mult_seq);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = LuzanEDoubleSparseMatrixMultSeqPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, LuzanEDoubleSparseMatrixMultSeqPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace luzan_e_double_sparse_matrix_mult_seq
