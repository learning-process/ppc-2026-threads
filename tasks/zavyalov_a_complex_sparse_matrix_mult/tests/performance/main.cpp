#include <gtest/gtest.h>

#include <cstddef>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

#include "util/include/perf_test_util.hpp"
#include "zavyalov_a_complex_sparse_matrix_mult/common/include/common.hpp"
#include "zavyalov_a_complex_sparse_matrix_mult/seq/include/ops_seq.hpp"

namespace zavyalov_a_compl_sparse_matr_mult {  // comm for ci1

class ZavyalovAComplexSparseMatrMultPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr size_t kCount = 15000;
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

    SparseMatrix expected(matr1 * matr2);

    return CompareSparse(expected, output_data);
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
