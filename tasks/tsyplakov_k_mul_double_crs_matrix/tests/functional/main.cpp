#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "tsyplakov_k_mul_double_crs_matrix/common/include/common.hpp"
#include "tsyplakov_k_mul_double_crs_matrix/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace tsyplakov_k_mul_double_crs_matrix {

class TsyplakovKRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    int case_num = std::get<0>(params);

    if (case_num == 1) {
      SparseMatrixCRS a(2, 2);
      a.values = {1.0, 1.0};
      a.col_index = {0, 1};
      a.row_ptr = {0, 1, 2};

      SparseMatrixCRS b(2, 2);
      b.values = {2.0, 3.0};
      b.col_index = {0, 1};
      b.row_ptr = {0, 1, 2};

      input_data_ = {a, b};
      expected_nnz_ = 2;
      expected_values_ = {2.0, 3.0};
    } else if (case_num == 2) {
      SparseMatrixCRS a(2, 2);
      a.values = {1.0, 2.0, 3.0};
      a.col_index = {0, 0, 1};
      a.row_ptr = {0, 1, 3};

      SparseMatrixCRS b(2, 2);
      b.values = {4.0, 5.0, 6.0};
      b.col_index = {0, 0, 1};
      b.row_ptr = {0, 1, 3};

      input_data_ = {a, b};
      expected_nnz_ = 3;
      expected_values_ = {4.0, 23.0, 18.0};
    } else if (case_num == 3) {
      SparseMatrixCRS a(2, 2);
      a.values = {1.0, 2.0};
      a.col_index = {0, 1};
      a.row_ptr = {0, 1, 2};

      SparseMatrixCRS b(2, 2);

      input_data_ = {a, b};
      expected_nnz_ = 0;
    } else if (case_num == 4) {
      SparseMatrixCRS a(3, 3);
      a.values = {1.0, 2.0, 3.0, 4.0};
      a.col_index = {0, 2, 1, 2};
      a.row_ptr = {0, 2, 3, 4};

      SparseMatrixCRS b(3, 3);
      b.values = {5.0, 6.0, 7.0, 8.0};
      b.col_index = {0, 1, 2, 0};
      b.row_ptr = {0, 1, 3, 4};

      input_data_ = {a, b};
      expected_nnz_ = 4;
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.rows != input_data_.a.rows || output_data.cols != input_data_.b.cols) {
      return false;
    }

    if (output_data.row_ptr.size() != static_cast<size_t>(output_data.rows + 1)) {
      return false;
    }

    if (output_data.values.size() != output_data.col_index.size()) {
      return false;
    }

    for (size_t i = 0; i < output_data.row_ptr.size() - 1; ++i) {
      if (output_data.row_ptr[i] > output_data.row_ptr[i + 1]) {
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
  int expected_nnz_ = 0;
  std::vector<double> expected_values_;
};

namespace {

TEST_P(TsyplakovKRunFuncTestsThreads, MatmulTest) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 4> kTestParam = {std::make_tuple(1, "identity"), std::make_tuple(2, "simple_2x2"),
                                            std::make_tuple(3, "zero_matrix"), std::make_tuple(4, "sparse_3x3")};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<TsyplakovKTestTaskSEQ, InType>(kTestParam, PPC_SETTINGS_tsyplakov_k_mul_double_crs_matrix));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);
const auto kPerfTestName = TsyplakovKRunFuncTestsThreads::PrintFuncTestName<TsyplakovKRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(MatrixMulTests, TsyplakovKRunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace tsyplakov_k_mul_double_crs_matrix
