#include <gtest/gtest.h>

#include <cstddef>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

#include "util/include/perf_test_util.hpp"
#include "zavyalov_a_complex_sparse_matr_mult/common/include/common.hpp"
#include "zavyalov_a_complex_sparse_matr_mult/seq/include/ops_seq.hpp"

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

    size_t rows_a = matr1.height;
    size_t cols_a_rows_b = matr1.width;
    size_t cols_b = matr2.width;

    std::vector<std::vector<Complex>> matr_a(rows_a, std::vector<Complex>(cols_a_rows_b, Complex(0.0, 0.0)));
    std::vector<std::vector<Complex>> matr_b(cols_a_rows_b, std::vector<Complex>(cols_b, Complex(0.0, 0.0)));

    for (size_t idx = 0; idx < matr1.Count(); ++idx) {
      size_t row = matr1.row_ind[idx];
      size_t col = matr1.col_ind[idx];
      Complex val = matr1.val[idx];
      if (row < rows_a && col < cols_a_rows_b) {
        matr_a[row][col] = val;
      }
    }

    for (size_t idx = 0; idx < matr2.Count(); ++idx) {
      size_t row = matr2.row_ind[idx];
      size_t col = matr2.col_ind[idx];
      Complex val = matr2.val[idx];
      if (row < cols_a_rows_b && col < cols_b) {
        matr_b[row][col] = val;
      }
    }

    std::vector<std::vector<Complex>> matr_c(rows_a, std::vector<Complex>(cols_b, Complex(0.0, 0.0)));

    for (size_t i = 0; i < rows_a; ++i) {
      for (size_t j = 0; j < cols_b; ++j) {
        for (size_t k = 0; k < cols_a_rows_b; ++k) {
          matr_c[i][j] += (matr_a[i][k] * matr_b[k][j]);
        }
      }
    }

    SparseMatrix expected(matr_c);
    if (expected.Count() != output_data.Count()) {
      return false;
    }

    std::map<std::pair<size_t, size_t>, Complex> output_map;
    for (size_t idx = 0; idx < output_data.Count(); ++idx) {
      output_map[{output_data.row_ind[idx], output_data.col_ind[idx]}] = output_data.val[idx];
    }

    for (size_t idx = 0; idx < expected.Count(); ++idx) {
      auto key = std::make_pair(expected.row_ind[idx], expected.col_ind[idx]);
      auto it = output_map.find(key);
      if (it == output_map.end()) {
        return false;
      }
      if (!(expected.val[idx] == it->second)) {
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
    PPC_SETTINGS_zavyalov_a_complex_sparse_matr_mult);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = ZavyalovAComplexSparseMatrMultPerfTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, ZavyalovAComplexSparseMatrMultPerfTest, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace zavyalov_a_compl_sparse_matr_mult
