#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "safaryan_a_sparse_matrix_mult_crs_seq/common/include/common.hpp"
#include "safaryan_a_sparse_matrix_mult_crs_seq/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace safaryan_a_sparse_matrix_mult_crs_seq {

namespace {

CRSMatrix CreateMatrix(size_t rows, size_t cols, const std::vector<double> &values,
                       const std::vector<size_t> &col_indices, const std::vector<size_t> &row_ptr) {
  CRSMatrix matrix;
  matrix.rows = rows;
  matrix.cols = cols;
  matrix.values = values;
  matrix.col_indices = col_indices;
  matrix.row_ptr = row_ptr;
  matrix.nnz = values.size();
  return matrix;
}

}  // namespace

class SafaryanARunFuncTestsSEQ : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<3>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    const int test_case = std::get<0>(params);

    switch (test_case) {
      case 1: {
        const CRSMatrix a = CreateMatrix(3, 3, {1.5, 4.2, 3.7, 2.8, 5.1}, {0, 2, 1, 0, 2}, {0, 2, 3, 5});
        const CRSMatrix b = CreateMatrix(3, 3, {1.2, 2.3, 3.4}, {0, 1, 2}, {0, 1, 2, 3});

        expected_output_ = CreateMatrix(3, 3, {1.8, 5.04, 8.51, 9.52, 17.34}, {0, 2, 1, 0, 2}, {0, 2, 3, 5});
        input_data_ = std::make_pair(a, b);
        break;
      }

      case 2: {
        const CRSMatrix a = CreateMatrix(2, 3, {1.2, 2.5, 3.7}, {0, 1, 2}, {0, 1, 2});
        const CRSMatrix b = CreateMatrix(3, 2, {1.1, 3.3, 2.2}, {0, 0, 1}, {0, 2, 3});

        expected_output_ = CreateMatrix(2, 2, {1.32, 12.21, 5.5}, {0, 1, 0}, {0, 2, 3});
        input_data_ = std::make_pair(a, b);
        break;
      }

      case 3: {
        const CRSMatrix a = CreateMatrix(2, 2, {2.5}, {0}, {0, 1, 1});
        const CRSMatrix b = CreateMatrix(2, 2, {3.7}, {1}, {0, 1, 1});

        expected_output_ = CreateMatrix(2, 2, {}, {}, {0, 0, 0});
        input_data_ = std::make_pair(a, b);
        break;
      }

      case 4: {
        const CRSMatrix a = CreateMatrix(3, 3, {1.1, 4.2, 7.3, 2.4, 5.5, 8.6, 3.7, 6.8, 9.9},
                                         {0, 1, 2, 0, 1, 2, 0, 1, 2}, {0, 3, 6, 9});
        const CRSMatrix b = CreateMatrix(3, 3, {1.0, 1.0, 1.0}, {0, 1, 2}, {0, 1, 2, 3});

        expected_output_ = a;
        input_data_ = std::make_pair(a, b);
        break;
      }

      case 5: {
        const CRSMatrix a = CreateMatrix(4, 4, {1.5, 2.5, 3.5, 4.5}, {0, 1, 2, 3}, {0, 1, 2, 3, 4});
        const CRSMatrix b = CreateMatrix(4, 4, {5.5, 6.5, 7.5, 8.5}, {0, 1, 2, 3}, {0, 1, 2, 3, 4});

        expected_output_ = CreateMatrix(4, 4, {8.25, 16.25, 26.25, 38.25}, {0, 1, 2, 3}, {0, 1, 2, 3, 4});
        input_data_ = std::make_pair(a, b);
        break;
      }

      default:
        break;
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.rows != expected_output_.rows || output_data.cols != expected_output_.cols) {
      return false;
    }

    if (output_data.values.size() != expected_output_.values.size() ||
        output_data.col_indices.size() != expected_output_.col_indices.size()) {
      return false;
    }

    if (output_data.row_ptr.size() != expected_output_.row_ptr.size()) {
      return false;
    }

    for (size_t i = 0; i < output_data.row_ptr.size(); ++i) {
      if (output_data.row_ptr[i] != expected_output_.row_ptr[i]) {
        return false;
      }
    }

    const double epsilon = 1e-10;
    for (size_t i = 0; i < output_data.values.size(); ++i) {
      if (std::abs(output_data.values[i] - expected_output_.values[i]) > epsilon) {
        return false;
      }
      if (output_data.col_indices[i] != expected_output_.col_indices[i]) {
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
  CRSMatrix expected_output_;
};

namespace {

TEST_P(SafaryanARunFuncTestsSEQ, SparseMatrixMultiplication) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 5> kTestParam = {
    std::make_tuple(1, 0, 0, "simple_3x3"), std::make_tuple(2, 0, 0, "rectangular_2x3_3x2"),
    std::make_tuple(3, 0, 0, "zero_result"), std::make_tuple(4, 0, 0, "identity_matrix"),
    std::make_tuple(5, 0, 0, "diagonal_matrices")};

const auto kTestTasksList =
    ppc::util::AddFuncTask<SafaryanATaskSEQ, InType>(kTestParam, PPC_SETTINGS_safaryan_a_sparse_matrix_mult_crs_seq);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = SafaryanARunFuncTestsSEQ::PrintFuncTestName<SafaryanARunFuncTestsSEQ>;

INSTANTIATE_TEST_SUITE_P(SparseMatrixMultFixedTests, SafaryanARunFuncTestsSEQ, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace safaryan_a_sparse_matrix_mult_crs_seq
