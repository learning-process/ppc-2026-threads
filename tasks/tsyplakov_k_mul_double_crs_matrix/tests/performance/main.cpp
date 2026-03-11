#include <gtest/gtest.h>

#include <vector>

#include "tsyplakov_k_mul_double_crs_matrix/common/include/common.hpp"
#include "tsyplakov_k_mul_double_crs_matrix/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace tsyplakov_k_mul_double_crs_matrix {

class TsyplakovKRunPerfTestsThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kSize_ = 100;

  InType input_data_{};

  void SetUp() override {
    SparseMatrixCRS a(kSize_, kSize_);
    SparseMatrixCRS b(kSize_, kSize_);

    std::vector<double> values_a;
    std::vector<double> values_b;
    std::vector<int> col_idx_a;
    std::vector<int> col_idx_b;
    std::vector<int> row_ptr_a(kSize_ + 1, 0);
    std::vector<int> row_ptr_b(kSize_ + 1, 0);

    int nnz_a = 0;
    int nnz_b = 0;
    for (int i = 0; i < kSize_; ++i) {
      row_ptr_a[i] = nnz_a;
      row_ptr_b[i] = nnz_b;

      values_a.push_back(1.0);
      col_idx_a.push_back(i);
      nnz_a++;

      values_b.push_back(1.0);
      col_idx_b.push_back(i);
      nnz_b++;
    }
    row_ptr_a[kSize_] = nnz_a;
    row_ptr_b[kSize_] = nnz_b;

    a.values = values_a;
    a.col_index = col_idx_a;
    a.row_ptr = row_ptr_a;

    b.values = values_b;
    b.col_index = col_idx_b;
    b.row_ptr = row_ptr_b;

    input_data_ = {.a = a, .b = b};
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return (output_data.rows == kSize_ && output_data.cols == kSize_ && !output_data.values.empty());
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(TsyplakovKRunPerfTestsThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, TsyplakovKTestTaskSEQ>(PPC_SETTINGS_tsyplakov_k_mul_double_crs_matrix);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = TsyplakovKRunPerfTestsThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, TsyplakovKRunPerfTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace tsyplakov_k_mul_double_crs_matrix
