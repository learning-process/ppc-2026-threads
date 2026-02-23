#include <gtest/gtest.h>

#include <algorithm>
#include <complex>
#include <random>
#include <tuple>
#include <vector>

#include "borunov_v_complex_ccs/common/include/common.hpp"
#include "borunov_v_complex_ccs/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace borunov_v_complex_ccs_seq {

SparseMatrix GenerateRandomSparseMatrix(int num_rows, int num_cols, double sparsity) {
  SparseMatrix mat;
  mat.num_rows = num_rows;
  mat.num_cols = num_cols;
  mat.col_ptrs.assign(num_cols + 1, 0);

  std::mt19937 gen(42);
  std::uniform_real_distribution<double> dist_val(-10.0, 10.0);
  std::uniform_real_distribution<double> dist_prob(0.0, 1.0);

  for (int j = 0; j < num_cols; ++j) {
    for (int i = 0; i < num_rows; ++i) {
      if (dist_prob(gen) < sparsity) {
        mat.values.push_back({dist_val(gen), dist_val(gen)});
        mat.row_indices.push_back(i);
      }
    }
    mat.col_ptrs[j + 1] = mat.values.size();
  }
  return mat;
}

SparseMatrix MultiplyDense(const SparseMatrix& A, const SparseMatrix& B) {
  SparseMatrix C;
  C.num_rows = A.num_rows;
  C.num_cols = B.num_cols;
  C.col_ptrs.assign(C.num_cols + 1, 0);

  std::vector<std::vector<std::complex<double>>> dense_A(A.num_rows, std::vector<std::complex<double>>(A.num_cols, {0.0, 0.0}));
  std::vector<std::vector<std::complex<double>>> dense_B(B.num_rows, std::vector<std::complex<double>>(B.num_cols, {0.0, 0.0}));
  std::vector<std::vector<std::complex<double>>> dense_C(A.num_rows, std::vector<std::complex<double>>(B.num_cols, {0.0, 0.0}));

  for (int j = 0; j < A.num_cols; ++j) {
    for (int idx = A.col_ptrs[j]; idx < A.col_ptrs[j + 1]; ++idx) {
      dense_A[A.row_indices[idx]][j] = A.values[idx];
    }
  }

  for (int j = 0; j < B.num_cols; ++j) {
    for (int idx = B.col_ptrs[j]; idx < B.col_ptrs[j + 1]; ++idx) {
      dense_B[B.row_indices[idx]][j] = B.values[idx];
    }
  }

  for (int i = 0; i < A.num_rows; ++i) {
    for (int j = 0; j < B.num_cols; ++j) {
      for (int k = 0; k < A.num_cols; ++k) {
        dense_C[i][j] += dense_A[i][k] * dense_B[k][j];
      }
    }
  }

  for (int j = 0; j < B.num_cols; ++j) {
    for (int i = 0; i < A.num_rows; ++i) {
      if (std::abs(dense_C[i][j]) > 1e-9) {
        C.values.push_back(dense_C[i][j]);
        C.row_indices.push_back(i);
      }
    }
    C.col_ptrs[j + 1] = C.values.size();
  }

  return C;
}

class BorunovVRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::to_string(std::get<1>(test_param)) + "_" + std::to_string(std::get<2>(test_param));
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int m = std::get<0>(params);
    int k = std::get<1>(params);
    int n = std::get<2>(params);

    SparseMatrix A = GenerateRandomSparseMatrix(m, k, 0.2);
    SparseMatrix B = GenerateRandomSparseMatrix(k, n, 0.2);

    input_data_ = {A, B};
    expected_output_ = {MultiplyDense(A, B)};
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (expected_output_.size() != output_data.size()) return false;
    const auto& expected = expected_output_[0];
    const auto& actual = output_data[0];

    if (expected.num_rows != actual.num_rows || expected.num_cols != actual.num_cols) return false;
    if (expected.col_ptrs != actual.col_ptrs) return false;
    if (expected.row_indices != actual.row_indices) return false;
    if (expected.values.size() != actual.values.size()) return false;

    for (size_t i = 0; i < expected.values.size(); ++i) {
      if (std::abs(expected.values[i] - actual.values[i]) > 1e-6) return false;
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_output_;
};

namespace {

TEST_P(BorunovVRunFuncTestsThreads, MatmulTests) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {std::make_tuple(10, 10, 10), std::make_tuple(20, 15, 25), std::make_tuple(5, 30, 5)};

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<BorunovVComplexCcsSEQ, InType>(kTestParam, PPC_SETTINGS_borunov_v_complex_ccs));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = BorunovVRunFuncTestsThreads::PrintFuncTestName<BorunovVRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(SparseMatrixTests, BorunovVRunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace borunov_v_complex_ccs_seq
