#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "kotelnikova_a_double_matr_mult/common/include/common.hpp"
#include "kotelnikova_a_double_matr_mult/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace kotelnikova_a_double_matr_mult {

namespace {
SparseMatrixCCS CreateMatrix(int rows, int cols, const std::vector<double> &values, const std::vector<int> &row_indices,
                             const std::vector<int> &col_ptrs) {
  SparseMatrixCCS matrix(rows, cols);
  matrix.values = values;
  matrix.row_indices = row_indices;
  matrix.col_ptrs = col_ptrs;
  return matrix;
}
}  // namespace

class KotelnikovaARunFuncTestsSEQ : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
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
        const SparseMatrixCCS a = CreateMatrix(3, 3, {1.5, 4.2, 3.7, 2.8, 5.1}, {0, 2, 1, 0, 2}, {0, 2, 3, 5});
        const SparseMatrixCCS b = CreateMatrix(3, 3, {1.2, 2.3, 3.4}, {0, 1, 2}, {0, 1, 2, 3});

        expected_output_ = CreateMatrix(3, 3, {1.8, 5.04, 8.51, 9.52, 17.34}, {0, 2, 1, 0, 2}, {0, 2, 3, 5});
        input_data_ = std::make_pair(a, b);
        break;
      }

      case 2: {
        const SparseMatrixCCS a = CreateMatrix(2, 3, {1.2, 2.5, 3.7}, {0, 0, 1}, {0, 1, 2, 3});
        const SparseMatrixCCS b = CreateMatrix(3, 2, {1.1, 3.3, 2.2}, {0, 2, 1}, {0, 2, 3});

        expected_output_ = CreateMatrix(2, 2, {1.32, 12.21, 5.5}, {0, 1, 0}, {0, 2, 3});
        input_data_ = std::make_pair(a, b);
        break;
      }

      case 3: {
        const SparseMatrixCCS a = CreateMatrix(2, 2, {2.5}, {0}, {0, 1, 1});
        const SparseMatrixCCS b = CreateMatrix(2, 2, {3.7}, {1}, {0, 1, 1});

        expected_output_ = CreateMatrix(2, 2, {}, {}, {0, 0, 0});
        input_data_ = std::make_pair(a, b);
        break;
      }

      case 4: {
        const SparseMatrixCCS a = CreateMatrix(3, 3, {1.1, 4.2, 7.3, 2.4, 5.5, 8.6, 3.7, 6.8, 9.9},
                                               {0, 1, 2, 0, 1, 2, 0, 1, 2}, {0, 3, 6, 9});
        const SparseMatrixCCS b = CreateMatrix(3, 3, {1.0, 1.0, 1.0}, {0, 1, 2}, {0, 1, 2, 3});

        expected_output_ = a;
        input_data_ = std::make_pair(a, b);
        break;
      }

      case 5: {
        const SparseMatrixCCS a = CreateMatrix(4, 4, {1.5, 2.5, 3.5, 4.5}, {0, 1, 2, 3}, {0, 1, 2, 3, 4});
        const SparseMatrixCCS b = CreateMatrix(4, 4, {5.5, 6.5, 7.5, 8.5}, {0, 1, 2, 3}, {0, 1, 2, 3, 4});

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
        output_data.row_indices.size() != expected_output_.row_indices.size()) {
      return false;
    }

    if (output_data.col_ptrs.size() != expected_output_.col_ptrs.size()) {
      return false;
    }

    for (size_t i = 0; i < output_data.col_ptrs.size(); ++i) {
      if (output_data.col_ptrs[i] != expected_output_.col_ptrs[i]) {
        return false;
      }
    }

    const double epsilon = 1e-10;
    for (size_t i = 0; i < output_data.values.size(); ++i) {
      if (std::abs(output_data.values[i] - expected_output_.values[i]) > epsilon) {
        return false;
      }
      if (output_data.row_indices[i] != expected_output_.row_indices[i]) {
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
  SparseMatrixCCS expected_output_;
};

namespace {

TEST_P(KotelnikovaARunFuncTestsSEQ, SparseMatrixMultiplication) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 5> kTestParam = {
    std::make_tuple(1, 0, 0, "simple_3x3"), std::make_tuple(2, 0, 0, "rectangular_2x3_3x2"),
    std::make_tuple(3, 0, 0, "zero_result"), std::make_tuple(4, 0, 0, "identity_matrix"),
    std::make_tuple(5, 0, 0, "diagonal_matrices")};

const auto kTestTasksList =
    ppc::util::AddFuncTask<KotelnikovaATaskSEQ, InType>(kTestParam, PPC_SETTINGS_kotelnikova_a_double_matr_mult);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = KotelnikovaARunFuncTestsSEQ::PrintFuncTestName<KotelnikovaARunFuncTestsSEQ>;

INSTANTIATE_TEST_SUITE_P(SparseMatrixMultFixedTests, KotelnikovaARunFuncTestsSEQ, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace kotelnikova_a_double_matr_mult
