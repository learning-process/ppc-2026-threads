#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "potashnik_m_matrix_mult_complex/common/include/common.hpp"
#include "potashnik_m_matrix_mult_complex/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace potashnik_m_matrix_mult_complex {

class PotashnikMMatrixMultComplexFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::to_string(std::get<1>(test_param)) + "_" +
           std::to_string(std::get<2>(test_param));
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());

    size_t rows_left = std::get<0>(params);
    size_t columns_left = std::get<1>(params);
    size_t rows_right = std::get<1>(params);
    size_t columns_right = std::get<2>(params);

    std::vector<std::vector<Complex>> matrix_left(rows_left);
    std::vector<std::vector<Complex>> matrix_right(rows_right);

    for (size_t i = 0; i < rows_left; ++i) {
      matrix_left[i].resize(columns_left);
      for (size_t j = 0; j < columns_left; ++j) {
        matrix_left[i][j] =
            Complex(static_cast<double>((i * i * rows_left) % 50U), static_cast<double>((j * j * columns_left) % 50U));
      }
    }

    for (size_t i = 0; i < rows_right; ++i) {
      matrix_right[i].resize(columns_right);
      for (size_t j = 0; j < columns_right; ++j) {
        matrix_right[i][j] = Complex(static_cast<double>((i * i * rows_right) % 50U),
                                     static_cast<double>((j * j * columns_right) % 50U));
      }
    }

    CCSMatrix matrix_left_CCS(matrix_left);
    CCSMatrix matrix_right_CCS(matrix_right);

    input_data_ = std::make_tuple(matrix_left_CCS, matrix_right_CCS);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    const CCSMatrix &matrix_left = std::get<0>(input_data_);
    const CCSMatrix &matrix_right = std::get<1>(input_data_);

    std::vector<Complex> val_left = matrix_left.val;
    std::vector<size_t> row_ind_left = matrix_left.row_ind;
    std::vector<size_t> col_ptr_left = matrix_left.col_ptr;
    size_t height_left = matrix_left.height;
    // size_t width_left = matrix_left.width;

    std::vector<Complex> val_right = matrix_right.val;
    std::vector<size_t> row_ind_right = matrix_right.row_ind;
    std::vector<size_t> col_ptr_right = matrix_right.col_ptr;
    // size_t height_right = matrix_right.height;
    size_t width_right = matrix_right.width;

    std::map<std::pair<size_t, size_t>, Complex> buffer;

    for (size_t i = 0; i < matrix_left.Count(); i++) {
      size_t row_left = row_ind_left[i];
      size_t col_left = col_ptr_left[i];
      Complex left_val = val_left[i];

      for (size_t j = 0; j < matrix_right.Count(); j++) {
        size_t row_right = row_ind_right[j];
        size_t col_right = col_ptr_right[j];
        Complex right_val = val_right[j];

        if (col_left == row_right) {
          buffer[{row_left, col_right}] += left_val * right_val;
        }
      }
    }

    CCSMatrix matrix_res;

    matrix_res.width = width_right;
    matrix_res.height = height_left;
    for (const auto &[key, value] : buffer) {
      matrix_res.val.push_back(value);
      matrix_res.row_ind.push_back(key.first);
      matrix_res.col_ptr.push_back(key.second);
    }

    return output_data.compare(matrix_res);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

namespace {

TEST_P(PotashnikMMatrixMultComplexFuncTests, MatmulFromPic) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 10> kTestParam = {
    std::make_tuple(4, 5, 6), std::make_tuple(6, 5, 4),   std::make_tuple(4, 6, 5), std::make_tuple(5, 6, 4),
    std::make_tuple(6, 6, 5), std::make_tuple(6, 5, 6),   std::make_tuple(5, 6, 6), std::make_tuple(8, 8, 8),
    std::make_tuple(9, 9, 9), std::make_tuple(10, 10, 10)};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<PotashnikMMatrixMultComplexSEQ, InType>(
    kTestParam, PPC_SETTINGS_potashnik_m_matrix_mult_complex));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName =
    PotashnikMMatrixMultComplexFuncTests::PrintFuncTestName<PotashnikMMatrixMultComplexFuncTests>;

INSTANTIATE_TEST_SUITE_P(FuncTests, PotashnikMMatrixMultComplexFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace potashnik_m_matrix_mult_complex
