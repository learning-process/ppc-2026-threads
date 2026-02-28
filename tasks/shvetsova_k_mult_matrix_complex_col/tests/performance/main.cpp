#include <gtest/gtest.h>

#include "shvetsova_k_mult_matrix_complex_col/common/include/common.hpp"
#include "shvetsova_k_mult_matrix_complex_col/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace shvetsova_k_mult_matrix_complex_col {

class ShvetsovaKRunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_{};

  void SetUp() override {
    const int size = 100;
    const int nnz_per_col = 5;

    MatrixCCS A;
    MatrixCCS B;

    A.rows = size;
    A.cols = size;
    A.col_ptr.resize(size + 1, 0);

    B.rows = size;
    B.cols = size;
    B.col_ptr.resize(size + 1, 0);

    // Матрица A
    for (int j = 0; j < size; ++j) {
      for (int k = 0; k < nnz_per_col; ++k) {
        int row = (j * nnz_per_col + k) % size;
        A.row_ind.push_back(row);
        A.values.emplace_back(1.0, 0.0);
      }
      A.col_ptr[j + 1] = static_cast<int>(A.values.size());
    }

    // Матрица B
    for (int j = 0; j < size; ++j) {
      for (int k = 0; k < nnz_per_col; ++k) {
        int row = (j * nnz_per_col + k) % size;
        B.row_ind.push_back(row);
        B.values.emplace_back(1.0, 0.0);
      }
      B.col_ptr[j + 1] = static_cast<int>(B.values.size());
    }

    input_data_ = std::make_tuple(A, B);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    const auto &matrix_C = output_data;
    return (matrix_C.cols > 0 && matrix_C.rows > 0) && static_cast<int>(matrix_C.col_ptr.size()) == matrix_C.cols + 1;
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
