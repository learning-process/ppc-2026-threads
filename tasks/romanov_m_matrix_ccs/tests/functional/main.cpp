#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "romanov_m_matrix_ccs/common/include/common.hpp"
#include "romanov_m_matrix_ccs/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace romanov_m_matrix_ccs {

class RomanovMRunFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int test_case = std::get<0>(params);

    MatrixCCS a;
    MatrixCCS b;
    if (test_case == 1) {
      a.rows_num = 2;
      a.cols_num = 2;
      a.nnz = 2;
      a.col_ptrs = {0, 1, 2};
      a.row_inds = {0, 1};
      a.vals = {2.0, 3.0};
      b = a;
      b.vals = {4.0, 5.0};
      expected_.rows_num = 2;
      expected_.cols_num = 2;
      expected_.nnz = 2;
      expected_.col_ptrs = {0, 1, 2};
      expected_.row_inds = {0, 1};
      expected_.vals = {8.0, 15.0};
    } else if (test_case == 2) {
      a.rows_num = 2;
      a.cols_num = 2;
      a.nnz = 2;
      a.col_ptrs = {0, 1, 2};
      a.row_inds = {0, 1};
      a.vals = {1.0, 1.0};
      b.rows_num = 2;
      b.cols_num = 2;
      b.nnz = 0;
      b.col_ptrs = {0, 0, 0};
      expected_.rows_num = 2;
      expected_.cols_num = 2;
      expected_.nnz = 0;
      expected_.col_ptrs = {0, 0, 0};
    } else if (test_case == 3) {
      a.rows_num = 1;
      a.cols_num = 3;
      a.nnz = 3;
      a.col_ptrs = {0, 1, 2, 3};
      a.row_inds = {0, 0, 0};
      a.vals = {1.0, 2.0, 3.0};
      b.rows_num = 3;
      b.cols_num = 1;
      b.nnz = 3;
      b.col_ptrs = {0, 3};
      b.row_inds = {0, 1, 2};
      b.vals = {4.0, 5.0, 6.0};
      expected_.rows_num = 1;
      expected_.cols_num = 1;
      expected_.nnz = 1;
      expected_.col_ptrs = {0, 1};
      expected_.row_inds = {0};
      expected_.vals = {32.0};
    }

    input_data_ = std::make_pair(a, b);
  }

  bool CheckTestOutputData(OutType &out) final {
    if (out.rows_num != expected_.rows_num || out.cols_num != expected_.cols_num) {
      return false;
    }
    if (out.vals.size() != expected_.vals.size()) {
      return false;
    }
    for (size_t i = 0; i < out.vals.size(); ++i) {
      if (std::abs(out.vals[i] - expected_.vals[i]) > 1e-9) {
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
  MatrixCCS matrix_a_, matrix_b_, expected_;
};

namespace {

TEST_P(RomanovMRunFuncTests, MatrixMultiplyTests) {
  ExecuteTest(GetParam());
}
const std::array<TestType, 3> kTestParam = {std::make_tuple(1, "test_2x2_basic"),
                                            std::make_tuple(2, "test_zero_matrix"),
                                            std::make_tuple(3, "test_rectangular_1x3_3x1")};
const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<RomanovMMatrixCCSSeq, InType>(kTestParam, PPC_SETTINGS_romanov_m_matrix_ccs));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = RomanovMRunFuncTests::PrintFuncTestName<RomanovMRunFuncTests>;

INSTANTIATE_TEST_SUITE_P(RomanovTests, RomanovMRunFuncTests, kGtestValues, kPerfTestName);

}  // namespace
}  // namespace romanov_m_matrix_ccs
