#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "agafonov_i_matrix_ccs_seq/common/include/common.hpp"
#include "agafonov_i_matrix_ccs_seq/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace agafonov_i_matrix_ccs_seq {

namespace {
CCSMatrix CreateMatrix(size_t rows, size_t cols, const std::vector<double> &values) {
  CCSMatrix matrix;
  matrix.rows_num = rows;
  matrix.cols_num = cols;
  matrix.col_ptrs.push_back(0);
  for (size_t j = 0; j < cols; ++j) {
    for (size_t i = 0; i < rows; ++i) {
      double val = values[(i * cols) + j];
      if (std::abs(val) > 1e-15) {
        matrix.vals.push_back(val);
        matrix.row_inds.push_back(static_cast<int>(i));
      }
    }
    matrix.col_ptrs.push_back(static_cast<int>(matrix.vals.size()));
  }
  matrix.nnz = matrix.vals.size();
  return matrix;
}
}  // namespace

class AgafonovIFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(
      const testing::TestParamInfo<ppc::util::BaseRunFuncTests<InType, OutType, TestType>::ParamType> &info) {
    return std::get<1>(std::get<static_cast<int>(ppc::util::GTestParamIndex::kTestParams)>(info.param));
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int test_id = std::get<0>(params);
    CCSMatrix matrix_a;
    CCSMatrix matrix_b;

    if (test_id == 0) {
      matrix_a = CreateMatrix(2, 2, {1, 0, 0, 2});
      matrix_b = CreateMatrix(2, 2, {3, 0, 0, 4});
      expected_ = CreateMatrix(2, 2, {3, 0, 0, 8});
    } else if (test_id == 1) {
      matrix_a = CreateMatrix(2, 2, {5, 6, 7, 8});
      matrix_b = CreateMatrix(2, 2, {1, 0, 0, 1});
      expected_ = matrix_a;
    } else if (test_id == 2) {
      matrix_a = CreateMatrix(2, 2, {1, 0, 1, 0});
      matrix_b = CreateMatrix(2, 2, {0, 0, 1, 1});
      expected_ = CreateMatrix(2, 2, {0, 0, 0, 0});
    } else if (test_id == 3) {
      matrix_a = CreateMatrix(1, 2, {1, 2});
      matrix_b = CreateMatrix(2, 1, {3, 4});
      expected_ = CreateMatrix(1, 1, {11});
    }

    input_data_ = std::make_pair(matrix_a, matrix_b);
  }

  bool CheckTestOutputData(OutType &out) final {
    if (out.rows_num != expected_.rows_num || out.cols_num != expected_.cols_num) {
      return false;
    }
    if (out.vals.size() != expected_.vals.size()) {
      return false;
    }
    if (out.col_ptrs != expected_.col_ptrs || out.row_inds != expected_.row_inds) {
      return false;
    }
    for (size_t i = 0; i < out.vals.size(); ++i) {
      if (std::abs(out.vals[i] - expected_.vals[i]) > 1e-12) {
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
  OutType expected_;
};

TEST_P(AgafonovIFuncTests, RunTests) {
  ExecuteTest(GetParam());
}

namespace {
const std::array<TestType, 4> kTestParams = {std::make_tuple(0, "Basic_2x2"), std::make_tuple(1, "Identity_Check"),
                                             std::make_tuple(2, "Zero_Result"),
                                             std::make_tuple(3, "Rectangular_Check")};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<AgafonovIMatrixCCSSeq, InType>(kTestParams, PPC_SETTINGS_agafonov_i_matrix_ccs_seq));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

INSTANTIATE_TEST_SUITE_P(AgafonovIndependentTests, AgafonovIFuncTests, kGtestValues,
                         AgafonovIFuncTests::PrintTestParam);
}  // namespace

}  // namespace agafonov_i_matrix_ccs_seq
