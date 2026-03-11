#include <gtest/gtest.h>

#include <cstddef>
#include <random>
#include <utility>
#include <vector>

#include "agafonov_i_matrix_ccs_seq/common/include/common.hpp"
#include "agafonov_i_matrix_ccs_seq/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace agafonov_i_matrix_ccs_seq {

namespace {
CCSMatrix CreateRandomCcs(size_t rows, size_t cols, double density) {
  CCSMatrix matrix;
  matrix.rows_num = rows;
  matrix.cols_num = cols;
  matrix.col_ptrs.push_back(0);

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> prob_dist(0.0, 1.0);
  std::uniform_real_distribution<double> val_dist(-10.0, 10.0);

  size_t nnz_count = 0;
  for (size_t j = 0; j < cols; ++j) {
    for (size_t i = 0; i < rows; ++i) {
      if (prob_dist(gen) < density) {
        matrix.vals.push_back(val_dist(gen));
        matrix.row_inds.push_back(static_cast<int>(i));
        nnz_count++;
      }
    }
    matrix.col_ptrs.push_back(static_cast<int>(matrix.vals.size()));
  }
  matrix.nnz = nnz_count;
  return matrix;
}
}  // namespace

class AgafonovMPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    size_t size = 1000;
    double density = 0.01;

    matrix_a_ = CreateRandomCcs(size, size, density);
    matrix_b_ = CreateRandomCcs(size, size, density);

    input_data_ = std::make_pair(matrix_a_, matrix_b_);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data.rows_num == matrix_a_.rows_num && output_data.cols_num == matrix_b_.cols_num;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  CCSMatrix matrix_a_;
  CCSMatrix matrix_b_;
  InType input_data_;
};

TEST_P(AgafonovMPerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {
const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, AgafonovIMatrixCCSSeq>(PPC_SETTINGS_agafonov_i_matrix_ccs_seq);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

INSTANTIATE_TEST_SUITE_P(MatrixPerfTests, AgafonovMPerfTest, kGtestValues, AgafonovMPerfTest::CustomPerfTestName);
}  // namespace

}  // namespace agafonov_i_matrix_ccs_seq
