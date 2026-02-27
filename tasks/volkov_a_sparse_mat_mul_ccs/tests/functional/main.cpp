#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"
#include "volkov_a_sparse_mat_mul_ccs/common/include/common.hpp"
#include "volkov_a_sparse_mat_mul_ccs/seq/include/ops_seq.hpp"

namespace volkov_a_sparse_mat_mul_ccs {

class VolkovAFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<0>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<int>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    std::string test_name = std::get<0>(params);

    SparseMatCCS mat_a;
    SparseMatCCS mat_b;

    if (test_name == "BasicMultiplication") {
      mat_a.rows_count = 3;
      mat_a.cols_count = 3;
      mat_a.non_zeros = 3;
      mat_a.col_ptrs = {0, 1, 2, 3};
      mat_a.row_indices = {0, 1, 2};
      mat_a.values = {1.0, 2.0, 3.0};

      mat_b.rows_count = 3;
      mat_b.cols_count = 2;
      mat_b.non_zeros = 4;
      mat_b.col_ptrs = {0, 2, 4};
      mat_b.row_indices = {0, 2, 0, 1};
      mat_b.values = {4.0, 7.0, 5.0, 6.0};

      expected_c_.rows_count = 3;
      expected_c_.cols_count = 2;
      expected_c_.non_zeros = 4;
      expected_c_.col_ptrs = {0, 2, 4};
      expected_c_.row_indices = {0, 2, 0, 1};
      expected_c_.values = {4.0, 21.0, 5.0, 12.0};

    } else if (test_name == "EmptyResultTest") {
      mat_a.rows_count = 2;
      mat_a.cols_count = 2;
      mat_a.non_zeros = 1;
      mat_a.col_ptrs = {0, 1, 1};
      mat_a.row_indices = {0};
      mat_a.values = {5.0};

      mat_b.rows_count = 2;
      mat_b.cols_count = 2;
      mat_b.non_zeros = 1;
      mat_b.col_ptrs = {0, 0, 1};
      mat_b.row_indices = {1};
      mat_b.values = {3.0};

      expected_c_.rows_count = 2;
      expected_c_.cols_count = 2;
      expected_c_.non_zeros = 0;
      expected_c_.col_ptrs = {0, 0, 0};
      expected_c_.row_indices = {};
      expected_c_.values = {};

    } else if (test_name == "VectorMultiplication") {
      mat_a.rows_count = 2;
      mat_a.cols_count = 2;
      mat_a.non_zeros = 4;
      mat_a.col_ptrs = {0, 2, 4};
      mat_a.row_indices = {0, 1, 0, 1};
      mat_a.values = {1.0, 3.0, 2.0, 4.0};

      mat_b.rows_count = 2;
      mat_b.cols_count = 1;
      mat_b.non_zeros = 2;
      mat_b.col_ptrs = {0, 2};
      mat_b.row_indices = {0, 1};
      mat_b.values = {5.0, 6.0};

      expected_c_.rows_count = 2;
      expected_c_.cols_count = 1;
      expected_c_.non_zeros = 2;
      expected_c_.col_ptrs = {0, 2};
      expected_c_.row_indices = {0, 1};
      expected_c_.values = {17.0, 39.0};

    } else if (test_name == "IdentityMatrixTest") {
      mat_a.rows_count = 2;
      mat_a.cols_count = 2;
      mat_a.non_zeros = 3;
      mat_a.col_ptrs = {0, 2, 3};
      mat_a.row_indices = {0, 1, 1};
      mat_a.values = {7.0, 8.0, 9.0};

      mat_b.rows_count = 2;
      mat_b.cols_count = 2;
      mat_b.non_zeros = 2;
      mat_b.col_ptrs = {0, 1, 2};
      mat_b.row_indices = {0, 1};
      mat_b.values = {1.0, 1.0};

      expected_c_ = mat_a;

    } else if (test_name == "NegativeValuesTest") {
      mat_a.rows_count = 2;
      mat_a.cols_count = 2;
      mat_a.non_zeros = 2;
      mat_a.col_ptrs = {0, 1, 2};
      mat_a.row_indices = {0, 1};
      mat_a.values = {-1.0, 2.0};

      mat_b.rows_count = 2;
      mat_b.cols_count = 2;
      mat_b.non_zeros = 3;
      mat_b.col_ptrs = {0, 1, 3};
      mat_b.row_indices = {0, 0, 1};
      mat_b.values = {3.0, -4.0, 5.0};

      expected_c_.rows_count = 2;
      expected_c_.cols_count = 2;
      expected_c_.non_zeros = 3;
      expected_c_.col_ptrs = {0, 1, 3};
      expected_c_.row_indices = {0, 0, 1};
      expected_c_.values = {-3.0, 4.0, 10.0};
    }

    input_data_ = std::make_tuple(mat_a, mat_b);
  }

  bool CheckTestOutputData(OutType &res) final {
    if (res.rows_count != expected_c_.rows_count) {
      return false;
    }
    if (res.cols_count != expected_c_.cols_count) {
      return false;
    }
    if (res.non_zeros != expected_c_.non_zeros) {
      return false;
    }

    if (res.col_ptrs != expected_c_.col_ptrs) {
      return false;
    }
    if (res.row_indices != expected_c_.row_indices) {
      return false;
    }

    for (size_t i = 0; i < expected_c_.values.size(); ++i) {
      if (std::abs(res.values[i] - expected_c_.values[i]) > 1e-9) {
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
  SparseMatCCS expected_c_;
};

namespace {

TEST_P(VolkovAFuncTests, RunIndependentTests) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 5> kTestParams = {
    std::make_tuple("BasicMultiplication", ""),
    std::make_tuple("EmptyResultTest", ""),
    std::make_tuple("VectorMultiplication", ""),
    std::make_tuple("IdentityMatrixTest", ""),
    std::make_tuple("NegativeValuesTest", "")
};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<VolkovASparseMatMulCcsSeq, InType>(
        kTestParams, PPC_SETTINGS_volkov_a_sparse_mat_mul_ccs));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);
const auto kPerfTestName = VolkovAFuncTests::PrintFuncTestName<VolkovAFuncTests>;

INSTANTIATE_TEST_SUITE_P(VolkovIndependentTests, VolkovAFuncTests, kGtestValues, kPerfTestName);

}  // namespace
}  // namespace volkov_a_sparse_mat_mul_ccs