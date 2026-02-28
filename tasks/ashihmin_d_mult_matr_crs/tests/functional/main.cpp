#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <string>
#include <tuple>
#include <vector>

#include "ashihmin_d_mult_matr_crs/common/include/common.hpp"
#include "ashihmin_d_mult_matr_crs/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace ashihmin_d_mult_matr_crs {

namespace {

CRSMatrix DenseToCRS(const DenseMatrix &dense_matrix) {
  CRSMatrix matrix_result;
  matrix_result.rows = static_cast<int>(dense_matrix.size());
  matrix_result.cols = dense_matrix.empty() ? 0 : static_cast<int>(dense_matrix[0].size());
  matrix_result.row_ptr.resize(matrix_result.rows + 1, 0);

  for (int row_index = 0; row_index < matrix_result.rows; ++row_index) {
    for (int col_index = 0; col_index < matrix_result.cols; ++col_index) {
      if (std::abs(dense_matrix[row_index][col_index]) > 1e-12) {
        matrix_result.values.push_back(dense_matrix[row_index][col_index]);
        matrix_result.col_index.push_back(col_index);
      }
    }
    matrix_result.row_ptr[row_index + 1] = static_cast<int>(matrix_result.values.size());
  }
  return matrix_result;
}

bool CompareCRS(const CRSMatrix &matrix_a, const CRSMatrix &matrix_b) {
  if (matrix_a.rows != matrix_b.rows || matrix_a.cols != matrix_b.cols) {
    return false;
  }
  if (matrix_a.row_ptr != matrix_b.row_ptr) {
    return false;
  }
  if (matrix_a.col_index != matrix_b.col_index) {
    return false;
  }
  if (matrix_a.values.size() != matrix_b.values.size()) {
    return false;
  }

  for (std::size_t index = 0; index < matrix_a.values.size(); ++index) {
    if (std::abs(matrix_a.values[index] - matrix_b.values[index]) > 1e-10) {
      return false;
    }
  }
  return true;
}

}  // namespace

class AshihminDMultMatrCrsFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &param) {
    return std::get<0>(param);
  }

 protected:
  void SetUp() override {
    TestType test_params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());

    const auto &dense_matrix_a = std::get<1>(test_params);
    const auto &dense_matrix_b = std::get<2>(test_params);
    const auto &dense_matrix_c = std::get<3>(test_params);

    input_data_ = std::make_tuple(DenseToCRS(dense_matrix_a), DenseToCRS(dense_matrix_b));
    expected_output_ = DenseToCRS(dense_matrix_c);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return CompareCRS(output_data, expected_output_);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  CRSMatrix expected_output_;
};

namespace {

TEST_P(AshihminDMultMatrCrsFuncTests, SpGemmCrsSeq) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 6> kTestParams = {
    std::make_tuple("Identity_2x2", DenseMatrix{{1, 0}, {0, 1}}, DenseMatrix{{5, 6}, {7, 8}},
                    DenseMatrix{{5, 6}, {7, 8}}),
    std::make_tuple("ZeroMatrix", DenseMatrix{{0, 0}, {0, 0}}, DenseMatrix{{1, 2}, {3, 4}},
                    DenseMatrix{{0, 0}, {0, 0}}),
    std::make_tuple("Simple_2x2", DenseMatrix{{1, 2}, {3, 4}}, DenseMatrix{{5, 6}, {7, 8}},
                    DenseMatrix{{19, 22}, {43, 50}}),
    std::make_tuple("Rectangular_2x3_3x2", DenseMatrix{{1, 0, 2}, {0, 3, 0}}, DenseMatrix{{0, 1}, {4, 0}, {5, 6}},
                    DenseMatrix{{10, 13}, {12, 0}}),
    std::make_tuple("SparseDiagonal", DenseMatrix{{1, 0, 0}, {0, 2, 0}, {0, 0, 3}},
                    DenseMatrix{{4, 0, 0}, {0, 5, 0}, {0, 0, 6}}, DenseMatrix{{4, 0, 0}, {0, 10, 0}, {0, 0, 18}}),
    std::make_tuple("SingleElement", DenseMatrix{{7}}, DenseMatrix{{8}}, DenseMatrix{{56}})};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<AshihminDMultMatrCrsSEQ, InType>(kTestParams, PPC_SETTINGS_ashihmin_d_mult_matr_crs));
const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = AshihminDMultMatrCrsFuncTests::PrintFuncTestName<AshihminDMultMatrCrsFuncTests>;
INSTANTIATE_TEST_SUITE_P(AshihminSparseCRSTests, AshihminDMultMatrCrsFuncTests, kGtestValues, kPerfTestName);

}  // namespace
}  // namespace ashihmin_d_mult_matr_crs
