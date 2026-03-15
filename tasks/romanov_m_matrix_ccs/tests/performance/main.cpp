#include <gtest/gtest.h>

#include <cstddef>
#include <random>
#include <utility>
#include <vector>

#include "romanov_m_matrix_ccs/common/include/common.hpp"
#include "romanov_m_matrix_ccs/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace romanov_m_matrix_ccs {

namespace {
MatrixCCS CreateRandomCCS(size_t rows, size_t cols, double density) {
  MatrixCCS matrix;
  matrix.rows_num = rows;
  matrix.cols_num = cols;
  matrix.col_ptrs.resize(cols + 1, 0);

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> prob_dist(0.0, 1.0);
  std::uniform_real_distribution<double> val_dist(-10.0, 10.0);

  size_t nnz_count = 0;
  for (size_t j = 0; j < cols; ++j) {
    matrix.col_ptrs[j] = nnz_count;
    for (size_t i = 0; i < rows; ++i) {
      if (prob_dist(gen) < density) {
        matrix.vals.push_back(val_dist(gen));
        matrix.row_inds.push_back(i);
        nnz_count++;
      }
    }
  }
  matrix.col_ptrs[cols] = nnz_count;
  matrix.nnz = nnz_count;
  return matrix;
}
}  // namespace

class RomanovMPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    size_t size = 1000;
    double density = 0.01;

    matrix_a_ = CreateRandomCCS(size, size, density);
    matrix_b_ = CreateRandomCCS(size, size, density);

    input_data_ = std::make_pair(matrix_a_, matrix_b_);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data.rows_num == matrix_a_.rows_num && output_data.cols_num == matrix_b_.cols_num;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  MatrixCCS matrix_a_;
  MatrixCCS matrix_b_;
  InType input_data_;
};

TEST_P(RomanovMPerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {
const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, RomanovMMatrixCCSSeq>(PPC_SETTINGS_romanov_m_matrix_ccs);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = RomanovMPerfTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(MatrixPerfTests, RomanovMPerfTest, kGtestValues, kPerfTestName);
}  // namespace

}  // namespace romanov_m_matrix_ccs
