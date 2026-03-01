#include <gtest/gtest.h>

#include <tuple>

#include "shvetsova_k_mult_matrix_complex_col/common/include/common.hpp"
#include "shvetsova_k_mult_matrix_complex_col/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace shvetsova_k_mult_matrix_complex_col {

class ShvetsovaKRunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;

  void SetUp() override {
    const int size = 100;
    const int nnz_per_col = 20;

    MatrixCCS matrix_a;
    MatrixCCS matrix_b;

    matrix_a.rows = size;
    matrix_a.cols = size;
    matrix_a.col_ptr.resize(size + 1, 0);

    matrix_b.rows = size;
    matrix_b.cols = size;
    matrix_b.col_ptr.resize(size + 1, 0);

    // Матрица matrix_a
    for (int j = 0; j < size; ++j) {
      for (int k = 0; k < nnz_per_col; ++k) {
        int row = (j * nnz_per_col + k) % size;
        matrix_a.row_ind.push_back(row);
        matrix_a.values.emplace_back(1.0, 0.0);
      }
      matrix_a.col_ptr[j + 1] = static_cast<int>(matrix_a.values.size());
    }

    // Матрица matrix_b
    for (int j = 0; j < size; ++j) {
      for (int k = 0; k < nnz_per_col; ++k) {
        int row = (j * nnz_per_col + k) % size;
        matrix_b.row_ind.push_back(row);
        matrix_b.values.emplace_back(1.0, 0.0);
      }
      matrix_b.col_ptr[j + 1] = static_cast<int>(matrix_b.values.size());
    }

    input_data_ = std::make_tuple(matrix_a, matrix_b);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return (output_data.cols > 0 && output_data.rows > 0);
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(ShvetsovaKRunPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, ShvetsovaKMultMatrixComplexSEQ>(
    PPC_SETTINGS_shvetsova_k_mult_matrix_complex_col);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = ShvetsovaKRunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, ShvetsovaKRunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace shvetsova_k_mult_matrix_complex_col
