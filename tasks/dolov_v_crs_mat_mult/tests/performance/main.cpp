#include <gtest/gtest.h>

#include <algorithm>
#include <tuple>
#include <vector>

#include "dolov_v_crs_mat_mult/common/include/common.hpp"
#include "dolov_v_crs_mat_mult/omp/include/ops_omp.hpp"
#include "dolov_v_crs_mat_mult/seq/include/ops_seq.hpp"
#include "dolov_v_crs_mat_mult/stl/include/ops_stl.hpp"
#include "dolov_v_crs_mat_mult/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace dolov_v_crs_mat_mult {

namespace {
SparseMatrix CreateBandMatrix(int n, int band_width) {
  SparseMatrix matrix;
  matrix.num_rows = n;
  matrix.num_cols = n;
  matrix.row_pointers.assign(n + 1, 0);

  for (int i = 0; i < n; ++i) {
    int start = std::max(0, i - band_width);
    int end = std::min(n - 1, i + band_width);
    for (int j = start; j <= end; ++j) {
      matrix.values.push_back(1.0);
      matrix.col_indices.push_back(j);
    }
    matrix.row_pointers[i + 1] = static_cast<int>(matrix.values.size());
  }
  return matrix;
}

}  // namespace

class DolovVCrsMatMultRunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    const int n = 2000;
    const int width = 20;

    SparseMatrix matrix_a = CreateBandMatrix(n, width);
    SparseMatrix matrix_b = CreateBandMatrix(n, width);

    input_data_ = {matrix_a, matrix_b};
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data.num_rows == input_data_[0].num_rows && !output_data.values.empty();
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

TEST_P(DolovVCrsMatMultRunPerfTestThreads, BandMatrixPerformance) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    std::tuple_cat(ppc::util::MakeAllPerfTasks<InType, DolovVCrsMatMultSeq>(PPC_SETTINGS_dolov_v_crs_mat_mult),
                   ppc::util::MakeAllPerfTasks<InType, DolovVCrsMatMultOmp>(PPC_SETTINGS_dolov_v_crs_mat_mult),
                   ppc::util::MakeAllPerfTasks<InType, DolovVCrsMatMultTbb>(PPC_SETTINGS_dolov_v_crs_mat_mult),
                   ppc::util::MakeAllPerfTasks<InType, DolovVCrsMatMultStl>(PPC_SETTINGS_dolov_v_crs_mat_mult));

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);
const auto kPerfTestName = DolovVCrsMatMultRunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(CRS_Band_Perf, DolovVCrsMatMultRunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace
}  // namespace dolov_v_crs_mat_mult
