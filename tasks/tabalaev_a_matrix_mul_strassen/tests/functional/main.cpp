#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>
#include <cmath>

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
    int r_a = std::get<0>(params);
    int c_a_r_b = std::get<1>(params);
    int c_b = std::get<2>(params);
    int up_to = std::get<3>(params);

    input_data_.a_rows = r_a;
    input_data_.a_cols_b_rows = c_a_r_b;
    input_data_.b_cols = c_b;
    input_data_.a.assign(r_a * c_a_r_b, 0.0);
    input_data_.b.assign(c_a_r_b * c_b, 0.0);

    for (int i = 0; i < r_a * c_a_r_b; ++i) {
      input_data_.a[i] = static_cast<double>(i % up_to);
    }
    for (int i = 0; i < c_a_r_b * c_b; ++i) {
      input_data_.b[i] = static_cast<double>(i % up_to) * 0.5;
    }

    expected_output_.assign(r_a * c_b, 0.0);
    for (int i = 0; i < r_a; ++i) {
      for (int k = 0; k < c_a_r_b; ++k) {
        double temp = input_data_.a[(i * c_a_r_b) + k];
        for (int j = 0; j < c_b; ++j) {
          expected_output_[(i * c_b) + j] += temp * input_data_.b[(k * c_b) + j];
        }
      }
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (output_data.size() != expected_output_.size()) {
      return false;
    }
    constexpr double kEps = 1e-9;
    for (size_t i = 0; i < output_data.size(); ++i) {
      if (std::fabs(output_data[i] - expected_output_[i]) > kEps) {
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
