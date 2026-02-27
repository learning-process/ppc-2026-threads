#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "tabalaev_a_matrix_mul_strassen/common/include/common.hpp"
#include "tabalaev_a_matrix_mul_strassen/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace tabalaev_a_matrix_mul_strassen {

class TabalaevAMatrixMulStrassenFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<4>(test_param) + "_" + std::to_string(std::get<0>(test_param)) + "x" +
           std::to_string(std::get<1>(test_param)) + "_" + std::to_string(std::get<1>(test_param)) + "x" +
           std::to_string(std::get<2>(test_param)) + "_Elems_Up_To_" + std::to_string(std::get<3>(test_param));
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    size_t rows_a = std::get<0>(params);
    size_t cols_a_rows_b = std::get<1>(params);
    size_t cols_b = std::get<2>(params);
    int up_to = std::get<3>(params);

    input_data_.a_rows = rows_a;
    input_data_.a_cols_b_rows = cols_a_rows_b;
    input_data_.b_cols = cols_b;
    input_data_.a.assign(rows_a * cols_a_rows_b, 0.0);
    input_data_.b.assign(cols_a_rows_b * cols_b, 0.0);

    for (size_t i = 0; i < rows_a * cols_a_rows_b; ++i) {
      input_data_.a[i] = static_cast<double>(i % up_to);
    }
    for (size_t i = 0; i < cols_a_rows_b * cols_b; ++i) {
      input_data_.b[i] = static_cast<double>(i % up_to) * 0.5;
    }

    expected_output_.assign(rows_a * cols_b, 0.0);
    for (size_t i = 0; i < rows_a; ++i) {
      for (size_t k = 0; k < cols_a_rows_b; ++k) {
        double temp = input_data_.a[(i * cols_a_rows_b) + k];
        for (size_t j = 0; j < cols_b; ++j) {
          expected_output_[(i * cols_b) + j] += temp * input_data_.b[(k * cols_b) + j];
        }
      }
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != expected_output_.size()) {
      return false;
    }
    constexpr double epsilon = 1e-9;
    for (size_t i = 0; i < output_data.size(); ++i) {
      if (std::fabs(output_data[i] - expected_output_[i]) > epsilon) {
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
  OutType expected_output_;
};

namespace {

TEST_P(TabalaevAMatrixMulStrassenFuncTests, MatmulFromPic) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 5> kTestParam = {
    std::make_tuple(2, 2, 2, 10, "Small_2x2"), std::make_tuple(4, 4, 4, 10, "PowerOfTwo_4x4"),
    std::make_tuple(3, 5, 3, 20, "NonSquare_Padded"), std::make_tuple(16, 16, 16, 100, "Medium_16x16"),
    std::make_tuple(64, 64, 64, 150, "Large_64x64")};

const auto kTestTasksList = ppc::util::AddFuncTask<TabalaevAMatrixMulStrassenSEQ, InType>(
    kTestParam, PPC_SETTINGS_tabalaev_a_matrix_mul_strassen);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = TabalaevAMatrixMulStrassenFuncTests::PrintFuncTestName<TabalaevAMatrixMulStrassenFuncTests>;

INSTANTIATE_TEST_SUITE_P(PicMatrixTests, TabalaevAMatrixMulStrassenFuncTests, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace tabalaev_a_matrix_mul_strassen
