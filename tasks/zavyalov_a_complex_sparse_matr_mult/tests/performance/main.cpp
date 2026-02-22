#include <gtest/gtest.h>

#include "util/include/perf_test_util.hpp"
#include "zavyalov_a_complex_sparse_matr_mult/common/include/common.hpp"
#include "zavyalov_a_complex_sparse_matr_mult/seq/include/ops_seq.hpp"

namespace zavyalov_a_compl_sparse_matr_mult {

class ZayalovAComplexSparseMatrMultPerfTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
  const size_t kCount_ = 1000;
  InType input_data_{};

  void SetUp() override {
    size_t n = kCount_;
    size_t m = kCount_;
    size_t k = kCount_;

    std::vector<std::vector<Complex>> matr_a(n, std::vector<Complex>(m, Complex(0.0, 0.0)));
    for (size_t i = 0; i < n; i++) {
      matr_a[i][(i * 43247u) % m] = Complex(43.0, 74.0);
      matr_a[i][(i * 73299u) % m] = Complex(i * 9.0, 7843.0);
    }

    std::vector<std::vector<Complex>> matr_b(m, std::vector<Complex>(k, Complex(0.0, 0.0)));
    for (size_t i = 0; i < m; i++) {
      matr_a[i][(i * 34627u) % k] = Complex(763.0, 743.0);
      matr_a[i][(i * 13337u) % k] = Complex(i * 953.0, 43215.0);
    }
    Sparse_matrix matr1(matr_a);
    Sparse_matrix matr2(matr_b);

    input_data_ = std::make_tuple(matr1, matr2);
  }

  bool CheckTestOutputData(OutType& output_data) final {
    Sparse_matrix& matr1 = std::get<0>(input_data_);
    Sparse_matrix& matr2 = std::get<1>(input_data_);
    size_t n = matr1.height;
    size_t m = matr1.width;
    size_t k = matr2.width;
    std::vector<std::vector<Complex>> matr_a(n);
    for (size_t i = 0; i < n; i++) {
      matr_a[i].resize(m);
    }
    std::vector<std::vector<Complex>> matr_b(m);
    for (size_t i = 0; i < m; i++) {
      matr_b[i].resize(k);
    }

    for (size_t i = 0; i < matr1.count(); i++) {
      size_t row = matr1.row_ind[i];
      size_t col = matr1.col_ind[i];
      Complex val = matr1.val[i];
      matr_a[row][col] = val;
    }

    for (size_t i = 0; i < matr2.count(); i++) {
      size_t row = matr2.row_ind[i];
      size_t col = matr2.col_ind[i];
      Complex val = matr2.val[i];
      matr_b[row][col] = val;
    }
    std::vector<std::vector<Complex>> matr_c(n, std::vector<Complex>(k, Complex(0, 0)));

    for (size_t i = 0; i < n; i++) {
      for (size_t j = 0; j < k; j++) {
        for (size_t p = 0; p < m; p++) {
          matr_c[i][j] += (matr_a[i][p] * matr_b[p][j]);
        }
      }
    }

    Sparse_matrix res(matr_c);
    if (res.count() != output_data.count()) {
      return false;
    }

    for (size_t i = 0; i < res.count(); i++) {
      bool ok = true;
      ok &= (res.row_ind[i] == output_data.row_ind[i]);
      ok &= (res.col_ind[i] == output_data.col_ind[i]);
      ok &= (res.val[i] == output_data.val[i]);
      if (!ok) {
        return false;
      }
    }

    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(ZayalovAComplexSparseMatrMultPerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, ZavyalovAComplSparseMatrMultSEQ>(
    PPC_SETTINGS_zavyalov_a_complex_sparse_matr_mult);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = ZayalovAComplexSparseMatrMultPerfTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, ZayalovAComplexSparseMatrMultPerfTest, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace zavyalov_a_compl_sparse_matr_mult
