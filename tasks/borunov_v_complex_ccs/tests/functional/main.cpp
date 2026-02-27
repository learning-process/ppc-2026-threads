#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <random>
#include <string>
#include <tuple>
#include <vector>

#include "borunov_v_complex_ccs/common/include/common.hpp"
#include "borunov_v_complex_ccs/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace borunov_v_complex_ccs_seq {

namespace {

using DenseMatrix = std::vector<std::vector<std::complex<double>>>;

SparseMatrix GenerateRandomSparseMatrix(int num_rows, int num_cols, double sparsity) {
  SparseMatrix mat;
  mat.num_rows = num_rows;
  mat.num_cols = num_cols;
  mat.col_ptrs.assign(num_cols + 1, 0);

  const int sparsity_seed = static_cast<int>(sparsity * 1000.0);
  std::seed_seq seed{num_rows, num_cols, sparsity_seed};
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> dist_val(-10.0, 10.0);
  std::uniform_real_distribution<double> dist_prob(0.0, 1.0);

  for (int j = 0; j < num_cols; ++j) {
    for (int i = 0; i < num_rows; ++i) {
      if (dist_prob(gen) < sparsity) {
        mat.values.emplace_back(dist_val(gen), dist_val(gen));
        mat.row_indices.push_back(i);
      }
    }
    mat.col_ptrs[j + 1] = static_cast<int>(mat.values.size());
  }
  return mat;
}

DenseMatrix BuildDense(const SparseMatrix &mat) {
  DenseMatrix dense(mat.num_rows, std::vector<std::complex<double>>(mat.num_cols, {0.0, 0.0}));
  for (int j = 0; j < mat.num_cols; ++j) {
    for (int idx = mat.col_ptrs[j]; idx < mat.col_ptrs[j + 1]; ++idx) {
      dense[mat.row_indices[idx]][j] = mat.values[idx];
    }
  }
  return dense;
}

DenseMatrix MultiplyDenseMatrices(const DenseMatrix &left, const DenseMatrix &right) {
  const int rows = static_cast<int>(left.size());
  const int inner = rows == 0 ? 0 : static_cast<int>(left[0].size());
  const int cols = right.empty() ? 0 : static_cast<int>(right[0].size());
  DenseMatrix result(rows, std::vector<std::complex<double>>(cols, {0.0, 0.0}));

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      for (int k = 0; k < inner; ++k) {
        result[i][j] += left[i][k] * right[k][j];
      }
    }
  }

  return result;
}

SparseMatrix BuildSparseFromDense(const DenseMatrix &dense) {
  SparseMatrix c;
  c.num_rows = static_cast<int>(dense.size());
  c.num_cols = dense.empty() ? 0 : static_cast<int>(dense[0].size());
  c.col_ptrs.assign(c.num_cols + 1, 0);

  for (int j = 0; j < c.num_cols; ++j) {
    for (int i = 0; i < c.num_rows; ++i) {
      if (std::abs(dense[i][j]) > 1e-9) {
        c.values.push_back(dense[i][j]);
        c.row_indices.push_back(i);
      }
    }
    c.col_ptrs[j + 1] = static_cast<int>(c.values.size());
  }

  return c;
}

SparseMatrix MultiplyDense(const SparseMatrix &a, const SparseMatrix &b) {
  DenseMatrix dense_a = BuildDense(a);
  DenseMatrix dense_b = BuildDense(b);
  DenseMatrix dense_c = MultiplyDenseMatrices(dense_a, dense_b);
  return BuildSparseFromDense(dense_c);
}

}  // namespace

class BorunovVRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::to_string(std::get<1>(test_param)) + "_" +
           std::to_string(std::get<2>(test_param));
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int m = std::get<0>(params);
    int k = std::get<1>(params);
    int n = std::get<2>(params);

    SparseMatrix a = GenerateRandomSparseMatrix(m, k, 0.2);
    SparseMatrix b = GenerateRandomSparseMatrix(k, n, 0.2);

    input_data_ = {a, b};
    expected_output_ = {MultiplyDense(a, b)};
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (expected_output_.size() != output_data.size()) {
      return false;
    }
    const auto &expected = expected_output_[0];
    const auto &actual = output_data[0];

    if (expected.num_rows != actual.num_rows || expected.num_cols != actual.num_cols) {
      return false;
    }
    if (expected.col_ptrs != actual.col_ptrs) {
      return false;
    }
    if (expected.row_indices != actual.row_indices) {
      return false;
    }
    if (expected.values.size() != actual.values.size()) {
      return false;
    }

    for (std::size_t i = 0; i < expected.values.size(); ++i) {
      if (std::abs(expected.values[i] - actual.values[i]) > 1e-6) {
        return false;
      }
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

const std::array<TestType, 3> kTestParam = {
    std::make_tuple(10, 10, 10),
    std::make_tuple(20, 15, 25),
    std::make_tuple(5, 30, 5),
};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<BorunovVComplexCcsSEQ, InType>(kTestParam, PPC_SETTINGS_borunov_v_complex_ccs));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = BorunovVRunFuncTestsThreads::PrintFuncTestName<BorunovVRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(SparseMatrixTests, BorunovVRunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace borunov_v_complex_ccs_seq
