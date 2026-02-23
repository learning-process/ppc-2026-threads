#include <gtest/gtest.h>

#include <algorithm>
#include <vector>

#include "borunov_v_complex_ccs/common/include/common.hpp"
#include "borunov_v_complex_ccs/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace borunov_v_complex_ccs_seq {

namespace {

SparseMatrix GenerateSparseMatrixPerf(int num_rows, int num_cols, int non_zeros_per_col) {
  SparseMatrix mat;
  mat.num_rows = num_rows;
  mat.num_cols = num_cols;
  mat.col_ptrs.assign(num_cols + 1, 0);

  int step = std::max(1, num_rows / non_zeros_per_col);

  for (int j = 0; j < num_cols; ++j) {
    for (int i = 0; i < num_rows; i += step) {
      mat.values.emplace_back(static_cast<double>(i + 1), static_cast<double>(j + 1));
      mat.row_indices.push_back(i);
    }
    mat.col_ptrs[j + 1] = static_cast<int>(mat.values.size());
  }
  return mat;
}

}  // namespace

class BorunovVRunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;

  void SetUp() override {
    int m = 200000;
    int k = 200000;
    int n = 200000;

    SparseMatrix a = GenerateSparseMatrixPerf(m, k, 20);
    SparseMatrix b = GenerateSparseMatrixPerf(k, n, 20);

    input_data_ = {a, b};
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data.size() == 1;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(BorunovVRunPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, BorunovVComplexCcsSEQ>(PPC_SETTINGS_borunov_v_complex_ccs);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = BorunovVRunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, BorunovVRunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace borunov_v_complex_ccs_seq
