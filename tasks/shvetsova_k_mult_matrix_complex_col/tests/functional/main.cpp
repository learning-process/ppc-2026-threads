#include <gtest/gtest.h>
#include <stb/stb_image.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "shvetsova_k_mult_matrix_complex_col/common/include/common.hpp"
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
        b.values = {{5, 0}, {7, 0}};

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

TEST_P(ShvetsovaKRunFuncTestsThreads, MatmulFromPic) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {std::make_tuple(0, "0"), std::make_tuple(1, "1"), std::make_tuple(2, "2")};

const auto kTestTasksList = std::tuple_cat(ppc::util::AddFuncTask<ShvetsovaKMultMatrixComplexSEQ, InType>(
    kTestParam, PPC_SETTINGS_shvetsova_k_mult_matrix_complex_col));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = ShvetsovaKRunFuncTestsThreads::PrintFuncTestName<ShvetsovaKRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(ShvetsovaKRunFuncTestsThreads, ShvetsovaKRunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace shvetsova_k_mult_matrix_complex_col
