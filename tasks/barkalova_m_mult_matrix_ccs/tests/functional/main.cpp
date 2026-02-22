#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "barkalova_m_mult_matrix_ccs/common/include/common.hpp"
#include "barkalova_m_mult_matrix_ccs/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace barkalova_m_mult_matrix_ccs {

class BarkalovaMatrixMultiplyFixedTest : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return "test_" + std::to_string(std::get<0>(test_param));
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
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
        matrix_a.values = {Complex(1.0, 2.0), Complex(3.0, 4.0)};
        matrix_a.nnz = static_cast<int>(matrix_a.values.size());

        matrix_b.rows = 2;
        matrix_b.cols = 2;
        matrix_b.col_ptrs = {0, 1, 2};
        matrix_b.row_indices = {0, 1};
        matrix_b.values = {Complex(2.0, -1.0), Complex(1.0, 1.0)};
        matrix_b.nnz = static_cast<int>(matrix_b.values.size());

        expected_result_.rows = 2;
        expected_result_.cols = 2;
        expected_result_.col_ptrs = {0, 1, 2};
        expected_result_.row_indices = {0, 1};
        expected_result_.values = {Complex(4.0, 3.0), Complex(-1.0, 7.0)};
        expected_result_.nnz = static_cast<int>(expected_result_.values.size());
        break;
      }

      case 2: {
        matrix_a.rows = 3;
        matrix_a.cols = 3;
        matrix_a.col_ptrs = {0, 1, 2, 3};
        matrix_a.row_indices = {0, 1, 2};
        matrix_a.values = {Complex(1.0, 0.0), Complex(1.0, 0.0), Complex(1.0, 0.0)};
        matrix_a.nnz = static_cast<int>(matrix_a.values.size());

        matrix_b.rows = 3;
        matrix_b.cols = 3;
        matrix_b.col_ptrs = {0, 1, 2, 3};
        matrix_b.row_indices = {0, 1, 2};
        matrix_b.values = {Complex(2.0, 3.0), Complex(4.0, -1.0), Complex(5.0, 2.0)};
        matrix_b.nnz = static_cast<int>(matrix_b.values.size());

        expected_result_.rows = 3;
        expected_result_.cols = 3;
        expected_result_.col_ptrs = {0, 1, 2, 3};
        expected_result_.row_indices = {0, 1, 2};
        expected_result_.values = {Complex(2.0, 3.0), Complex(4.0, -1.0), Complex(5.0, 2.0)};
        expected_result_.nnz = static_cast<int>(expected_result_.values.size());
        break;
      }

      case 3: {
        matrix_a.rows = 2;
        matrix_a.cols = 2;
        matrix_a.col_ptrs = {0, 1, 2};
        matrix_a.row_indices = {0, 1};
        matrix_a.values = {Complex(1.0, 0.0), Complex(2.0, 0.0)};
        matrix_a.nnz = static_cast<int>(matrix_a.values.size());

        matrix_b.rows = 2;
        matrix_b.cols = 2;
        matrix_b.col_ptrs = {0, 1, 2};
        matrix_b.row_indices = {0, 1};
        matrix_b.values = {Complex(3.0, 0.0), Complex(4.0, 0.0)};
        matrix_b.nnz = static_cast<int>(matrix_b.values.size());

        expected_result_.rows = 2;
        expected_result_.cols = 2;
        expected_result_.col_ptrs = {0, 1, 2};
        expected_result_.row_indices = {0, 1};
        expected_result_.values = {Complex(3.0, 0.0), Complex(8.0, 0.0)};
        expected_result_.nnz = 2;
        break;
      }

      case 4: {
        matrix_a.rows = 2;
        matrix_a.cols = 2;
        matrix_a.col_ptrs = {0, 1, 2};
        matrix_a.row_indices = {0, 1};
        matrix_a.values = {Complex(1.0, 1.0), Complex(2.0, -1.0)};
        matrix_a.nnz = 2;

        matrix_b.rows = 2;
        matrix_b.cols = 3;
        matrix_b.col_ptrs = {0, 2, 4, 6};
        matrix_b.row_indices = {0, 1, 0, 1, 0, 1};
        matrix_b.values = {Complex(2.0, 0.0), Complex(5.0, 0.0), Complex(3.0, 0.0),
                           Complex(6.0, 0.0), Complex(4.0, 0.0), Complex(7.0, 0.0)};
        matrix_b.nnz = 6;

        expected_result_.rows = 2;
        expected_result_.cols = 3;
        expected_result_.col_ptrs = {0, 2, 4, 6};
        expected_result_.row_indices = {0, 1, 0, 1, 0, 1};
        expected_result_.values = {Complex(2.0, 2.0),   Complex(10.0, -5.0), Complex(3.0, 3.0),
                                   Complex(12.0, -6.0), Complex(4.0, 4.0),   Complex(14.0, -7.0)};
        expected_result_.nnz = 6;
        break;
      }

      case 5: {
        matrix_a.rows = 3;
        matrix_a.cols = 3;
        matrix_a.col_ptrs = {0, 2, 3, 4};
        matrix_a.row_indices = {0, 2, 1, 2};
        matrix_a.values = {Complex(1.0, 2.0), Complex(3.0, 1.0), Complex(2.0, -3.0), Complex(4.0, -2.0)};
        matrix_a.nnz = 4;

        matrix_b.rows = 3;
        matrix_b.cols = 3;
        matrix_b.col_ptrs = {0, 2, 3, 4};
        matrix_b.row_indices = {0, 2, 1, 2};
        matrix_b.values = {Complex(5.0, -1.0), Complex(7.0, 3.0), Complex(6.0, 2.0), Complex(8.0, -4.0)};
        matrix_b.nnz = 4;

        expected_result_.rows = 3;
        expected_result_.cols = 3;
        expected_result_.col_ptrs = {0, 2, 3, 4};
        expected_result_.row_indices = {0, 2, 1, 2};
        expected_result_.values = {Complex(7.0, 9.0), Complex(50.0, 0.0), Complex(18.0, -14.0), Complex(24.0, -32.0)};
        expected_result_.nnz = 4;
        break;
      }
      case 6: {
        matrix_a.rows = 3;
        matrix_a.cols = 2;
        matrix_a.col_ptrs = {0, 2, 4};
        matrix_a.row_indices = {0, 1, 0, 2};
        matrix_a.values = {Complex(1.0, 1.0), Complex(3.0, 3.0), Complex(2.0, -2.0), Complex(4.0, -4.0)};
        matrix_a.nnz = 4;

        matrix_b.rows = 2;
        matrix_b.cols = 4;
        matrix_b.col_ptrs = {0, 2, 4, 6, 8};
        matrix_b.row_indices = {0, 1, 0, 1, 0, 1, 0, 1};
        matrix_b.values = {Complex(1.0, 0.0), Complex(5.0, 0.0), Complex(2.0, 0.0), Complex(6.0, 0.0),
                           Complex(3.0, 0.0), Complex(7.0, 0.0), Complex(4.0, 0.0), Complex(8.0, 0.0)};
        matrix_b.nnz = 8;

        expected_result_.rows = 3;
        expected_result_.cols = 4;
        expected_result_.col_ptrs = {0, 3, 6, 9, 12};
        expected_result_.row_indices = {0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2};
        expected_result_.values = {Complex(11.0, -9.0),  Complex(3.0, 3.0),   Complex(20.0, -20.0),
                                   Complex(14.0, -10.0), Complex(6.0, 6.0),   Complex(24.0, -24.0),
                                   Complex(17.0, -11.0), Complex(9.0, 9.0),   Complex(28.0, -28.0),
                                   Complex(20.0, -12.0), Complex(12.0, 12.0), Complex(32.0, -32.0)};
        expected_result_.nnz = 12;
        break;
      }

      default:
        throw std::runtime_error("Unknown test case");
    }

    input_data_ = std::make_pair(matrix_a, matrix_b);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    const double eps = 1e-10;

    if (output_data.rows != expected_result_.rows || output_data.cols != expected_result_.cols) {
      return false;
    }

    if (output_data.nnz != expected_result_.nnz) {
      return false;
    }

    if (output_data.col_ptrs.size() != expected_result_.col_ptrs.size()) {
      return false;
    }

    for (size_t i = 0; i < output_data.col_ptrs.size(); ++i) {
      if (output_data.col_ptrs[i] != expected_result_.col_ptrs[i]) {
        return false;
      }
    }

    if (output_data.row_indices.size() != expected_result_.row_indices.size()) {
      return false;
    }

    for (size_t i = 0; i < output_data.row_indices.size(); ++i) {
      if (output_data.row_indices[i] != expected_result_.row_indices[i]) {
        return false;
      }
    }

    for (size_t i = 0; i < output_data.values.size(); ++i) {
      if (std::abs(output_data.values[i].real() - expected_result_.values[i].real()) > eps ||
          std::abs(output_data.values[i].imag() - expected_result_.values[i].imag()) > eps) {
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
};

namespace {

TEST_P(BarkalovaMatrixMultiplyFixedTest, MatrixMultiplyFixedTest) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 6> kFixedTestParams = {std::make_tuple(1, ""), std::make_tuple(2, ""),
                                                  std::make_tuple(3, ""), std::make_tuple(4, ""),
                                                  std::make_tuple(5, ""), std::make_tuple(6, "")};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<BarkalovaMMultMatrixCcsSEQ, InType>(
    kFixedTestParams, PPC_SETTINGS_barkalova_m_mult_matrix_ccs));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kTestName = BarkalovaMatrixMultiplyFixedTest::PrintFuncTestName<BarkalovaMatrixMultiplyFixedTest>;

INSTANTIATE_TEST_SUITE_P(FixedMatrixTests, BarkalovaMatrixMultiplyFixedTest, kGtestValues, kTestName);

}  // namespace
}  // namespace barkalova_m_mult_matrix_ccs
