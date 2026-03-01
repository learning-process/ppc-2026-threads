#include <gtest/gtest.h>
#include <algorithm>
#include <vector>
#include <compare>

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
    SparseMatrixCCS a, b;
    int n = kCount_;

    a.rows = a.cols = n;
    b.rows = b.cols = n;

    a.col_ptr.resize(n + 1, 0);
    b.col_ptr.resize(n + 1, 0);

    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<double> dist(-10.0, 10.0);

    for (int j = 0; j < n; j++) {
      std::vector<int> rows_a, rows_b;

      while (static_cast<int>(rows_a.size()) < std::min(kNonZeroPerCol_, n)) {
        int r = static_cast<int>(rng() % static_cast<std::mt19937::result_type>(n));
        if (std::find(rows_a.begin(), rows_a.end(), r) == rows_a.end()) {
          rows_a.push_back(r);
        }
      }
      while ((int)rows_b.size() < std::min(kNonZeroPerCol_, n)) {
        int r = rng() % n;
        if (std::find(rows_b.begin(), rows_b.end(), r) == rows_b.end()) {
          rows_b.push_back(r);
        }
      }

      std::sort(rows_a.begin(), rows_a.end());
      std::sort(rows_b.begin(), rows_b.end());

      for (int r : rows_a) {
        a.row_ind.push_back(r);
        a.values.emplace_back(dist(rng), dist(rng));
      }
      a.col_ptr[j + 1] = a.col_ptr[j] + rows_a.size();

      for (int r : rows_b) {
        b.row_ind.push_back(r);
        b.values.emplace_back(dist(rng), dist(rng));
      }
      b.col_ptr[j + 1] = b.col_ptr[j] + rows_b.size();
    }

    input_data_ = {a, b};
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
