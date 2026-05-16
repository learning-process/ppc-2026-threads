#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "kapanova_s_sparse_matrix_mult_ccs_seq/common/include/common.hpp"
#include "kapanova_s_sparse_matrix_mult_ccs_seq/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace kapanova_s_sparse_matrix_mult_ccs_seq {

class KapanovaSMatrixMultiplyTest : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return "test_" + std::to_string(std::get<0>(test_param));
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    params_ = params;
    int test_case = std::get<0>(params);

    CCSMatrix matrix_a;
    CCSMatrix matrix_b;
    expected_result_ = CCSMatrix();

    switch (test_case) {
      case 1: {
        matrix_a.rows = 2;
        matrix_a.cols = 2;
        matrix_a.col_ptrs = {0, 1, 2};
        matrix_a.row_indices = {0, 1};
        matrix_a.values = {2.0, 3.0};
        matrix_a.nnz = 2;

        matrix_b.rows = 2;
        matrix_b.cols = 2;
        matrix_b.col_ptrs = {0, 1, 2};
        matrix_b.row_indices = {0, 1};
        matrix_b.values = {4.0, 5.0};
        matrix_b.nnz = 2;

        expected_result_.rows = 2;
        expected_result_.cols = 2;
        expected_result_.col_ptrs = {0, 1, 2};
        expected_result_.row_indices = {0, 1};
        expected_result_.values = {8.0, 15.0};
        expected_result_.nnz = 2;
        break;
      }

      case 2: {
        matrix_a.rows = 3;
        matrix_a.cols = 3;
        matrix_a.col_ptrs = {0, 1, 3, 4};
        matrix_a.row_indices = {0, 1, 2, 2};
        matrix_a.values = {1.0, 2.0, 3.0, 4.0};
        matrix_a.nnz = 4;

        matrix_b.rows = 3;
        matrix_b.cols = 3;
        matrix_b.col_ptrs = {0, 2, 3, 5};
        matrix_b.row_indices = {0, 2, 1, 1, 2};
        matrix_b.values = {5.0, 6.0, 7.0, 8.0, 9.0};
        matrix_b.nnz = 5;

        expected_result_.rows = 3;
        expected_result_.cols = 3;
        expected_result_.col_ptrs = {0, 2, 4, 6};
        expected_result_.row_indices = {0, 2, 1, 2, 1, 2};
        expected_result_.values = {5.0, 24.0, 14.0, 21.0, 16.0, 60.0};
        expected_result_.nnz = 6;
        break;
      }

      case 3: {
        matrix_a.rows = 3;
        matrix_a.cols = 3;
        matrix_a.col_ptrs = {0, 2, 3, 4};
        matrix_a.row_indices = {0, 2, 1, 2};
        matrix_a.values = {1.0, 3.0, 2.0, 4.0};
        matrix_a.nnz = 4;

        matrix_b.rows = 3;
        matrix_b.cols = 3;
        matrix_b.col_ptrs = {0, 1, 3, 4};
        matrix_b.row_indices = {1, 0, 2, 2};
        matrix_b.values = {5.0, 6.0, 7.0, 8.0};
        matrix_b.nnz = 4;

        expected_result_.rows = 3;
        expected_result_.cols = 3;
        expected_result_.col_ptrs = {0, 1, 3, 4};
        expected_result_.row_indices = {1, 0, 2, 2};
        expected_result_.values = {10.0, 6.0, 46.0, 32.0};
        expected_result_.nnz = 4;
        break;
      }

      case 4: {
        matrix_a.rows = 2;
        matrix_a.cols = 3;
        matrix_a.col_ptrs = {0, 2, 3, 5};
        matrix_a.row_indices = {0, 1, 1, 0, 1};
        matrix_a.values = {1.0, 4.0, 2.0, 3.0, 5.0};
        matrix_a.nnz = 5;

        matrix_b.rows = 3;
        matrix_b.cols = 2;
        matrix_b.col_ptrs = {0, 2, 5};
        matrix_b.row_indices = {0, 2, 1, 2, 1};
        matrix_b.values = {6.0, 8.0, 7.0, 9.0, 10.0};
        matrix_b.nnz = 5;

        expected_result_.rows = 2;
        expected_result_.cols = 2;
        expected_result_.col_ptrs = {0, 2, 4};
        expected_result_.row_indices = {0, 1, 0, 1};
        expected_result_.values = {30.0, 64.0, 27.0, 79.0};
        expected_result_.nnz = 4;
        break;
      }

      case 5: {
        matrix_a.rows = 2;
        matrix_a.cols = 2;
        matrix_a.col_ptrs = {0, 0, 0};
        matrix_a.row_indices = {};
        matrix_a.values = {};
        matrix_a.nnz = 0;

        matrix_b.rows = 2;
        matrix_b.cols = 2;
        matrix_b.col_ptrs = {0, 1, 2};
        matrix_b.row_indices = {0, 1};
        matrix_b.values = {1.0, 2.0};
        matrix_b.nnz = 2;

        expected_result_.rows = 2;
        expected_result_.cols = 2;
        expected_result_.col_ptrs = {0, 0, 0};
        expected_result_.row_indices = {};
        expected_result_.values = {};
        expected_result_.nnz = 0;
        break;
      }

      default:
        throw std::runtime_error("Unknown test case");
    }

    input_data_ = std::make_pair(matrix_a, matrix_b);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    const double eps = 1e-10;

    if (output_data.rows != expected_result_.rows) {
      return false;
    }
    if (output_data.cols != expected_result_.cols) {
      return false;
    }
    if (output_data.nnz != expected_result_.nnz) {
      return false;
    }
    if (output_data.col_ptrs != expected_result_.col_ptrs) {
      return false;
    }
    if (output_data.row_indices != expected_result_.row_indices) {
      return false;
    }

    for (size_t i = 0; i < output_data.values.size(); ++i) {
      if (std::abs(output_data.values[i] - expected_result_.values[i]) > eps) {
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
  CCSMatrix expected_result_;
  TestType params_;
};

namespace {

TEST_P(KapanovaSMatrixMultiplyTest, MatrixMultiplyFixedTest) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 5> kFixedTestParams = {std::make_tuple(1, ""), std::make_tuple(2, ""),
                                                  std::make_tuple(3, ""), std::make_tuple(4, ""),
                                                  std::make_tuple(5, "")};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<KapanovaSSparseMatrixMultCCSSeq, InType>(
    kFixedTestParams, PPC_SETTINGS_kapanova_s_sparse_matrix_mult_ccs_seq));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kTestName = KapanovaSMatrixMultiplyTest::PrintFuncTestName<KapanovaSMatrixMultiplyTest>;

INSTANTIATE_TEST_SUITE_P(FixedMatrixTests, KapanovaSMatrixMultiplyTest, kGtestValues, kTestName);

}  // namespace

}  // namespace kapanova_s_sparse_matrix_mult_ccs_seq
