#include <gtest/gtest.h>

#include <cstddef>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

#include "util/include/perf_test_util.hpp"
#include "zavyalov_a_complex_sparse_matrix_mult/common/include/common.hpp"
#include "zavyalov_a_complex_sparse_matrix_mult/seq/include/ops_seq.hpp"

namespace zavyalov_a_compl_sparse_matr_mult {

class ZavyalovAComplexSparseMatrMultPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr size_t kCount = 1000;
  InType input_data_;

  void SetUp() override {
    size_t rows_a = kCount;
    size_t cols_a_rows_b = kCount;
    size_t cols_b = kCount;

    std::vector<std::vector<Complex>> matr_a(rows_a, std::vector<Complex>(cols_a_rows_b, Complex(0.0, 0.0)));
    for (size_t i = 0; i < rows_a; ++i) {
      matr_a[i][(i * 43247U) % cols_a_rows_b] = Complex(43.0, 74.0);
      matr_a[i][(i * 73299U) % cols_a_rows_b] = Complex(static_cast<double>(i) * 9.0, 7843.0);
    }

    std::vector<std::vector<Complex>> matr_b(cols_a_rows_b, std::vector<Complex>(cols_b, Complex(0.0, 0.0)));
    for (size_t i = 0; i < cols_a_rows_b; ++i) {
      matr_b[i][(i * 34627U) % cols_b] = Complex(763.0, 743.0);
      matr_b[i][(i * 13337U) % cols_b] = Complex(static_cast<double>(i) * 953.0, 43215.0);
    }

    SparseMatrix matr1(matr_a);
    SparseMatrix matr2(matr_b);

    input_data_ = std::make_tuple(matr1, matr2);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    const SparseMatrix &matr1 = std::get<0>(input_data_);
    const SparseMatrix &matr2 = std::get<1>(input_data_);

    auto matr_a = RestoreDense(matr1, matr1.height, matr1.width);
    auto matr_b = RestoreDense(matr2, matr2.height, matr2.width);

    auto matr_c = MultiplyDense(matr_a, matr_b);

    SparseMatrix expected(matr_c);

    return CompareSparse(expected, output_data);
  }

  std::vector<std::vector<Complex>> static RestoreDense(const SparseMatrix &matr, size_t rows, size_t cols) {
    std::vector<std::vector<Complex>> result(rows, std::vector<Complex>(cols, Complex(0.0, 0.0)));

    for (size_t idx = 0; idx < matr.Count(); ++idx) {
      size_t row = matr.row_ind[idx];
      size_t col = matr.col_ind[idx];
      if (row < rows && col < cols) {
        result[row][col] = matr.val[idx];
      }
    }

    return result;
  }
  std::vector<std::vector<Complex>> static MultiplyDense(const std::vector<std::vector<Complex>> &a,
                                                         const std::vector<std::vector<Complex>> &b) {
    size_t rows = a.size();
    size_t inner = b.size();
    size_t cols = b[0].size();

    std::vector<std::vector<Complex>> result(rows, std::vector<Complex>(cols, Complex(0.0, 0.0)));

    for (size_t i = 0; i < rows; ++i) {
      for (size_t j = 0; j < cols; ++j) {
        for (size_t k = 0; k < inner; ++k) {
          result[i][j] += a[i][k] * b[k][j];
        }
      }
    }

    return result;
  }
  bool static CompareSparse(const SparseMatrix &expected, const SparseMatrix &output) {
    if (expected.Count() != output.Count()) {
      return false;
    }

    std::map<std::pair<size_t, size_t>, Complex> output_map;

    for (size_t idx = 0; idx < output.Count(); ++idx) {
      output_map[{output.row_ind[idx], output.col_ind[idx]}] = output.val[idx];
    }

    for (size_t idx = 0; idx < expected.Count(); ++idx) {
      auto key = std::make_pair(expected.row_ind[idx], expected.col_ind[idx]);

      auto it = output_map.find(key);
      if (it == output_map.end() || !(expected.val[idx] == it->second)) {
        return false;
      }
    }

    return true;
  }
  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(ZavyalovAComplexSparseMatrMultPerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, ZavyalovAComplSparseMatrMultSEQ>(
    PPC_SETTINGS_zavyalov_a_complex_sparse_matrix_mult);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = ZavyalovAComplexSparseMatrMultPerfTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, ZavyalovAComplexSparseMatrMultPerfTest, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace zavyalov_a_compl_sparse_matr_mult
