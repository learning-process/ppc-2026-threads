#include <gtest/gtest.h>

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

#include "maslova_u_mult_matr_crs/common/include/common.hpp"
#include "maslova_u_mult_matr_crs/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace maslova_u_mult_matr_crs {

class MaslovaUMultMatrRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return "test_case_" + std::to_string(std::get<0>(test_param));
  }

 protected:
  void SetUp() override {
    int test_case =
        std::get<0>(std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam()));
    CRSMatrix A, B;
    expected_result_ = CRSMatrix();

    switch (test_case) {
      case 1: {
        A.rows = 2;
        A.cols = 2;
        A.row_ptr = {0, 2, 4};
        A.col_ind = {0, 1, 0, 1};
        A.values = {1.0, 1.0, 1.0, 1.0};
        B.rows = 2;
        B.cols = 2;
        B.row_ptr = {0, 2, 4};
        B.col_ind = {0, 1, 0, 1};
        B.values = {1.0, 1.0, 1.0, 1.0};

        expected_result_.rows = 2;
        expected_result_.cols = 2;
        expected_result_.row_ptr = {0, 2, 4};
        expected_result_.col_ind = {0, 1, 0, 1};
        expected_result_.values = {2.0, 2.0, 2.0, 2.0};
        break;
      }

      case 2: {
        A.rows = 2;
        A.cols = 2;
        A.row_ptr = {0, 1, 3};
        A.col_ind = {0, 0, 1};
        A.values = {1.0, 1.0, 1.0};
        B.rows = 2;
        B.cols = 2;
        B.row_ptr = {0, 1, 3};
        B.col_ind = {0, 0, 1};
        B.values = {1.0, 1.0, 1.0};

        expected_result_.rows = 2;
        expected_result_.cols = 2;
        expected_result_.row_ptr = {0, 1, 3};
        expected_result_.col_ind = {0, 0, 1};
        expected_result_.values = {1.0, 2.0, 1.0};
        break;
      }

      case 3: {
        A.rows = 1;
        A.cols = 1;
        A.row_ptr = {0, 1};
        A.col_ind = {0};
        A.values = {5.0};
        B.rows = 1;
        B.cols = 2;
        B.row_ptr = {0, 2};
        B.col_ind = {0, 1};
        B.values = {2.0, 3.0};

        expected_result_.rows = 1;
        expected_result_.cols = 2;
        expected_result_.row_ptr = {0, 2};
        expected_result_.col_ind = {0, 1};
        expected_result_.values = {10.0, 15.0};
        break;
      }

      case 4: {
        A.rows = 2;
        A.cols = 2;
        A.row_ptr = {0, 1, 1};
        A.col_ind = {1};
        A.values = {5.0};
        B.rows = 2;
        B.cols = 2;
        B.row_ptr = {0, 0, 1};
        B.col_ind = {0};
        B.values = {2.0};

        expected_result_.rows = 2;
        expected_result_.cols = 2;
        expected_result_.row_ptr = {0, 1, 1};
        expected_result_.col_ind = {0};
        expected_result_.values = {10.0};
        break;
      }

      case 5: {
        A.rows = 1;
        A.cols = 3;
        A.row_ptr = {0, 3};
        A.col_ind = {0, 1, 2};
        A.values = {1.0, 2.0, 3.0};
        B.rows = 3;
        B.cols = 1;
        B.row_ptr = {0, 1, 2, 3};
        B.col_ind = {0, 0, 0};
        B.values = {4.0, 5.0, 6.0};

        expected_result_.rows = 1;
        expected_result_.cols = 1;
        expected_result_.row_ptr = {0, 1};
        expected_result_.col_ind = {0};
        expected_result_.values = {32.0};
        break;
      }

      default:
        throw std::runtime_error("Invalid Test Case");
    }
    input_data_ = std::make_tuple(A, B);
  }

  bool CheckTestOutputData(OutType &out) final {
    if (std::tie(out.rows, out.cols, out.row_ptr, out.col_ind) !=
        std::tie(expected_result_.rows, expected_result_.cols, expected_result_.row_ptr, expected_result_.col_ind)) {
      return false;
    }

    if (out.values.size() != expected_result_.values.size()) {
      return false;
    }

    for (size_t i = 0; i < out.values.size(); ++i) {
      if (std::abs(out.values[i] - expected_result_.values[i]) > 1e-9) {
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
  OutType expected_result_;
};

namespace {

TEST_P(MaslovaUMultMatrRunFuncTestsThreads, MultMatrCRS) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 5> kTestParams = {std::make_tuple(1, "FullDense2x2"), std::make_tuple(2, "LowerTriangular"),
                                             std::make_tuple(3, "OuterProduct"), std::make_tuple(4, "CrossShift"),
                                             std::make_tuple(5, "VectorInnerProduct")};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<MaslovaUMultMatrSEQ, InType>(kTestParams, PPC_SETTINGS_maslova_u_mult_matr_crs));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = MaslovaUMultMatrRunFuncTestsThreads::PrintFuncTestName<MaslovaUMultMatrRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(MultMatrCRS, MaslovaUMultMatrRunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace maslova_u_mult_matr_crs
