#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "dilshodov_a_spmm_double_css/common/include/common.hpp"
#include "dilshodov_a_spmm_double_css/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace dilshodov_a_spmm_double_css {

namespace {
constexpr double kEps = 1e-10;

SparseMatrixCCS DenseToCCS(const DenseMatrix &dense) {
  SparseMatrixCCS matrix;
  matrix.rows_count = static_cast<int>(dense.size());
  matrix.cols_count = dense.empty() ? 0 : static_cast<int>(dense[0].size());
  matrix.col_ptrs.assign(static_cast<size_t>(matrix.cols_count) + 1, 0);

  for (int col = 0; col < matrix.cols_count; ++col) {
    for (int row = 0; row < matrix.rows_count; ++row) {
      if (std::abs(dense[row][col]) > kEps) {
        matrix.row_indices.push_back(row);
        matrix.values.push_back(dense[row][col]);
      }
    }
    matrix.col_ptrs[col + 1] = static_cast<int>(matrix.values.size());
  }

  matrix.non_zeros = static_cast<int>(matrix.values.size());
  return matrix;
}

DenseMatrix DenseMul(const DenseMatrix &lhs, const DenseMatrix &rhs) {
  const int rows = static_cast<int>(lhs.size());
  const int common = lhs.empty() ? 0 : static_cast<int>(lhs[0].size());
  const int cols = rhs.empty() ? 0 : static_cast<int>(rhs[0].size());

  DenseMatrix result(static_cast<size_t>(rows), std::vector<double>(static_cast<size_t>(cols), 0.0));
  for (int row = 0; row < rows; ++row) {
    for (int pivot = 0; pivot < common; ++pivot) {
      const double lhs_value = lhs[row][pivot];
      if (std::abs(lhs_value) <= kEps) {
        continue;
      }
      for (int col = 0; col < cols; ++col) {
        result[row][col] += lhs_value * rhs[pivot][col];
      }
    }
  }
  return result;
}

bool CompareCCS(const SparseMatrixCCS &lhs, const SparseMatrixCCS &rhs) {
  if (lhs.rows_count != rhs.rows_count || lhs.cols_count != rhs.cols_count || lhs.non_zeros != rhs.non_zeros) {
    return false;
  }
  if (lhs.col_ptrs != rhs.col_ptrs || lhs.row_indices != rhs.row_indices || lhs.values.size() != rhs.values.size()) {
    return false;
  }

  for (size_t i = 0; i < lhs.values.size(); ++i) {
    if (std::abs(lhs.values[i] - rhs.values[i]) > 1e-8) {
      return false;
    }
  }
  return true;
}

}  // namespace

class DilshodovASpmmDoubleCssExampleFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &param) {
    return std::get<0>(param);
  }

 protected:
  void SetUp() override {
    const TestType test_param = std::get<static_cast<size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    const auto &dense_a = std::get<1>(test_param);
    const auto &dense_b = std::get<2>(test_param);

    input_data_ = std::make_tuple(DenseToCCS(dense_a), DenseToCCS(dense_b));
    expected_output_ = DenseToCCS(DenseMul(dense_a, dense_b));
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return CompareCCS(output_data, expected_output_);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  SparseMatrixCCS expected_output_;
};

namespace {

TEST_P(DilshodovASpmmDoubleCssExampleFuncTests, RunExampleLikeTests) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParams = {
    std::make_tuple("TwoByTwoBasic", DenseMatrix{{1.0, 0.0}, {2.0, 3.0}}, DenseMatrix{{4.0, 5.0}, {0.0, 6.0}}),
    std::make_tuple("ThreeByThreeSparse", DenseMatrix{{0.0, 2.0, 0.0}, {1.0, 0.0, 3.0}, {0.0, 4.0, 0.0}},
                    DenseMatrix{{7.0, 0.0, 1.0}, {0.0, 8.0, 0.0}, {2.0, 0.0, 9.0}}),
    std::make_tuple("RectangularCheck", DenseMatrix{{1.0, -1.0, 0.0}, {0.0, 2.0, 4.0}},
                    DenseMatrix{{3.0, 0.0}, {0.0, 5.0}, {6.0, 7.0}})};

const auto kTestTasksList =
    ppc::util::AddFuncTask<DilshodovASpmmDoubleCssSeq, InType>(kTestParams, PPC_SETTINGS_dilshodov_a_spmm_double_css);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName =
    DilshodovASpmmDoubleCssExampleFuncTests::PrintFuncTestName<DilshodovASpmmDoubleCssExampleFuncTests>;

INSTANTIATE_TEST_SUITE_P(DilshodovSpmmDoubleCssExampleFunc, DilshodovASpmmDoubleCssExampleFuncTests, kGtestValues,
                         kPerfTestName);

}  // namespace

}  // namespace dilshodov_a_spmm_double_css
