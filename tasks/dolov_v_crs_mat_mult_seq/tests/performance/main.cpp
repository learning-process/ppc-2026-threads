#include <gtest/gtest.h>
#include <vector>
#include <numeric>

#include "dolov_v_crs_mat_mult_seq/common/include/common.hpp"
#include "dolov_v_crs_mat_mult_seq/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace dolov_v_crs_mat_mult_seq {

namespace {
SparseMatrix CreateBandMatrix(int n, int band_width) {
  SparseMatrix matrix;
  matrix.num_rows = n;
  matrix.num_cols = n;
  matrix.row_pointers.assign(n + 1, 0);

  for (int i = 0; i < n; ++i) {
    for (int j = std::max(0, i - band_width); j <= std::min(n - 1, i + band_width); ++j) {
      matrix.values.push_back(1.0);
      matrix.col_indices.push_back(j);
    }
    matrix.row_pointers[i + 1] = static_cast<int>(matrix.values.size());
  }
  return matrix;
}

}  // namespace

class DolovVCrsMatMultSeqRunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    const int n = 1000;
    const int width = 10;

    SparseMatrix a = CreateBandMatrix(n, width);
    SparseMatrix b = CreateBandMatrix(n, width);
    
    input_data_ = {a, b};
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data.num_rows == 1000 && !output_data.values.empty();
  }

  InType GetTestInputData() final { return input_data_; }

 private:
  InType input_data_;
};

TEST_P(DolovVCrsMatMultSeqRunPerfTestThreads, BandMatrixPerformance) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, DolovVCrsMatMultSeq>(PPC_SETTINGS_dolov_v_crs_mat_mult_seq);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = DolovVCrsMatMultSeqRunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(Sequential_Band_Perf, DolovVCrsMatMultSeqRunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace
}  // namespace dolov_v_crs_mat_mult_seq