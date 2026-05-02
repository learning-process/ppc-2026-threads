#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "../../common/include/common.hpp"
#include "../../omp/include/ops_omp.hpp"
#include "../../stl/include/ops_stl.hpp"
#include "../../tbb/include/ops_tbb.hpp"
#include "omp.h"
#include "shvetsova_k_mult_matrix_complex_col/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace shvetsova_k_mult_matrix_complex_col {

class ShvetsovaKRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());

    int test_id = std::get<0>(params);

    MatrixCCS a;
    MatrixCCS b;
    MatrixCCS c_expected;

    switch (test_id) {
      case 0: {
        a.rows = 2;
        a.cols = 2;
        a.col_ptr = {0, 1, 2};
        a.row_ind = {0, 1};
        a.values = {{1, 0}, {1, 0}};

        // B
        b.rows = 2;
        b.cols = 2;
        b.col_ptr = {0, 1, 2};
        b.row_ind = {0, 1};
        b.values = {{5, 2}, {7, 1}};

        // C = A*B = B
        c_expected = b;
        break;
      }

      case 1: {
        // A = diag(2,3)
        a.rows = 2;
        a.cols = 2;
        a.col_ptr = {0, 1, 2};
        a.row_ind = {0, 1};
        a.values = {{2, 0}, {3, 0}};

        // B = diag(4,5)
        b.rows = 2;
        b.cols = 2;
        b.col_ptr = {0, 1, 2};
        b.row_ind = {0, 1};
        b.values = {{4, 0}, {5, 0}};

        // C = diag(8,15)
        c_expected.rows = 2;
        c_expected.cols = 2;
        c_expected.col_ptr = {0, 1, 2};
        c_expected.row_ind = {0, 1};
        c_expected.values = {{8, 0}, {15, 0}};
        break;
      }

      case 2: {
        a.rows = 2;
        a.cols = 2;
        a.col_ptr = {0, 0, 0};

        b.rows = 2;
        b.cols = 2;
        b.col_ptr = {0, 0, 0};

        c_expected = a;
        break;
      }
      case 3: {
        // A = [[1+i, 0], [0, 2+2i]]
        a.rows = 2;
        a.cols = 2;
        a.col_ptr = {0, 1, 2};
        a.row_ind = {0, 1};
        a.values = {{1, 1}, {2, 2}};

        // B = [[0+1i, 1+0i], [3+4i, 0-1i]]
        b.rows = 2;
        b.cols = 2;
        b.col_ptr = {0, 2, 4};
        b.row_ind = {0, 1, 0, 1};
        b.values = {{0, 1}, {3, 4}, {1, 0}, {0, -1}};

        // C = A*B
        // C[0,0] = (1+i)*(0+1i) = -1+1i
        // C[1,0] = (2+2i)*(3+4i) = -2+14i
        // C[0,1] = (1+i)*(1+0i) = 1+1i
        // C[1,1] = (2+2i)*(0-1i) = 2-2i
        c_expected.rows = 2;
        c_expected.cols = 2;
        c_expected.col_ptr = {0, 2, 4};
        c_expected.row_ind = {0, 1, 0, 1};
        c_expected.values = {{-1, 1}, {-2, 14}, {1, 1}, {2, -2}};
        break;
      }

      case 4: {
        // A 1x2, B 2x1 -> результат 1x1
        a.rows = 1;
        a.cols = 2;
        a.col_ptr = {0, 1, 2};
        a.row_ind = {0, 0};
        a.values = {{1, 1}, {2, -1}};

        b.rows = 2;
        b.cols = 1;
        b.col_ptr = {0, 2};
        b.row_ind = {0, 1};
        b.values = {{0, 1}, {1, 1}};

        // C = A*B
        // C[0,0] = (1+i)*(0+1i) + (2-1i)*(1+1i) = (-1+1i) + (3+1i) = (2+2i)
        c_expected.rows = 1;
        c_expected.cols = 1;
        c_expected.col_ptr = {0, 1};
        c_expected.row_ind = {0};
        c_expected.values = {{2, 2}};
        break;
      }
      default:
        throw std::runtime_error("Unknown test_id");
    }

    input_data_ = {a, b};
    expected_data_ = {c_expected};
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return output_data == expected_data_;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_data_;
};

namespace {

TEST_P(ShvetsovaKRunFuncTestsThreads, MultMatrixComplex) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 5> kTestParam = {std::make_tuple(0, "0"), std::make_tuple(1, "1"), std::make_tuple(2, "2"),
                                            std::make_tuple(3, "3"), std::make_tuple(4, "4")};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<ShvetsovaKMultMatrixComplexSEQ, InType>(
                                               kTestParam, PPC_SETTINGS_shvetsova_k_mult_matrix_complex_col),
                                           ppc::util::AddFuncTask<ShvetsovaKMultMatrixComplexOMP, InType>(
                                               kTestParam, PPC_SETTINGS_shvetsova_k_mult_matrix_complex_col),
                                           ppc::util::AddFuncTask<ShvetsovaKMultMatrixComplexTBB, InType>(
                                               kTestParam, PPC_SETTINGS_shvetsova_k_mult_matrix_complex_col),
                                           ppc::util::AddFuncTask<ShvetsovaKMultMatrixComplexSTL, InType>(
                                               kTestParam, PPC_SETTINGS_shvetsova_k_mult_matrix_complex_col));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = ShvetsovaKRunFuncTestsThreads::PrintFuncTestName<ShvetsovaKRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(ShvetsovaKRunFuncTestsThreads, ShvetsovaKRunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace shvetsova_k_mult_matrix_complex_col
