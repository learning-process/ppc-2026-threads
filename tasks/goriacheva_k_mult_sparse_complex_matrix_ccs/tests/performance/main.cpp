#include <gtest/gtest.h>

#include <random>

#include "goriacheva_k_mult_sparse_complex_matrix_ccs/common/include/common.hpp"
#include "goriacheva_k_mult_sparse_complex_matrix_ccs/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace goriacheva_k_mult_sparse_complex_matrix_ccs {

class GoriachevaKMultSparseComplexMatrixCcsPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const int kCount_ = 200;
  const int kNonZeroPerCol_ = 5;
  InType input_data_{};

  void SetUp() override {
    SparseMatrixCCS A, B;
    int n = kCount_;

    A.rows = A.cols = n;
    B.rows = B.cols = n;

    A.col_ptr.resize(n + 1, 0);
    B.col_ptr.resize(n + 1, 0);

    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(-10.0, 10.0);

    for (int j = 0; j < n; j++) {
      std::vector<int> rows_A, rows_B;

      while ((int)rows_A.size() < std::min(kNonZeroPerCol_, n)) {
        int r = rng() % n;
        if (std::find(rows_A.begin(), rows_A.end(), r) == rows_A.end()) {
          rows_A.push_back(r);
        }
      }
      while ((int)rows_B.size() < std::min(kNonZeroPerCol_, n)) {
        int r = rng() % n;
        if (std::find(rows_B.begin(), rows_B.end(), r) == rows_B.end()) {
          rows_B.push_back(r);
        }
      }

      std::sort(rows_A.begin(), rows_A.end());
      std::sort(rows_B.begin(), rows_B.end());

      for (int r : rows_A) {
        A.row_ind.push_back(r);
        A.values.emplace_back(dist(rng), dist(rng));
      }
      A.col_ptr[j + 1] = A.col_ptr[j] + rows_A.size();

      for (int r : rows_B) {
        B.row_ind.push_back(r);
        B.values.emplace_back(dist(rng), dist(rng));
      }
      B.col_ptr[j + 1] = B.col_ptr[j] + rows_B.size();
    }

    input_data_ = {A, B};
  }

  bool CheckTestOutputData(OutType & /*output_data*/) final {
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(GoriachevaKMultSparseComplexMatrixCcsPerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, GoriachevaKMultSparseComplexMatrixCcsSEQ>(
    PPC_SETTINGS_goriacheva_k_mult_sparse_complex_matrix_ccs);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = GoriachevaKMultSparseComplexMatrixCcsPerfTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, GoriachevaKMultSparseComplexMatrixCcsPerfTest, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace goriacheva_k_mult_sparse_complex_matrix_ccs
