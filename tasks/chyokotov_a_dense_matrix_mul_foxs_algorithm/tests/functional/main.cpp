#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "chyokotov_a_dense_matrix_mul_foxs_algorithm/common/include/common.hpp"
#include "chyokotov_a_dense_matrix_mul_foxs_algorithm/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace chyokotov_a_dense_matrix_mul_foxs_algorithm {

class ChyokotovARunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    std::vector<double> a = std::get<0>(test_param).first;
    if (!a.empty()) {
      int n = static_cast<int>(std::sqrt(static_cast<double>(a.size())));
      return "size_of_output_matrix_" + std::to_string(n) + "x" + std::to_string(n);
    }
    return "empty_matrix";
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
    expected_output_ = std::get<1>(params);
    std::vector<double> a = input_data_.first;
    std::vector<double> b = input_data_.second;
    int n = static_cast<int>(std::sqrt(static_cast<double>(a.size())));

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        double sum = 0.0;
        for (int k = 0; k < n; k++) {
          sum += a[(i * n) + k] * b[(k * n) + j];
        }
        expected_output_[(i * n) + j] = sum;
      }
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    return (expected_output_ == output_data);
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  OutType expected_output_;
};

namespace {

TEST_P(ChyokotovARunFuncTestsThreads, MatrixMultiplicate) {
  ExecuteTest(GetParam());
}

const std::vector<double> kEmptymatrix = {};
const std::vector<double> kMatrixa1x1 = {1.0};
const std::vector<double> kMatrixb1x1 = {2.0};
const std::vector<double> kMatrixc1x1 = {2.0};

const std::vector<double> kMatrixa2x2 = {1.0, 2.0, 1.5, 2.5};
const std::vector<double> kMatrixb2x2 = {2.0, 4.0, 1.5, 3.5};
const std::vector<double> kMatrixc2x2 = {5.0, 11.0, 6.75, 14.75};

const std::vector<double> kMatrixa4x4 = {1.0, 2.0, 3.0, 4.0, 4.0, 3.0, 2.0, 1.0,
                                         5.0, 6.0, 7.0, 8.0, 8.0, 7.0, 6.0, 5.0};
const std::vector<double> kMatrixb4x4 = {1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                                         0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0};

const std::array<TestType, 4> kTestParam = {
    std::make_tuple(std::make_pair(kEmptymatrix, kEmptymatrix), kEmptymatrix),
    std::make_tuple(std::make_pair(kMatrixa1x1, kMatrixb1x1), kMatrixc1x1),
    std::make_tuple(std::make_pair(kMatrixa2x2, kMatrixb2x2), kMatrixc2x2),
    std::make_tuple(std::make_pair(kMatrixa4x4, kMatrixb4x4), kMatrixa4x4),
};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<ChyokotovADenseMatMulFoxAlgorithmSEQ, InType>(kTestParam, PPC_SETTINGS_example_threads));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = ChyokotovARunFuncTestsThreads::PrintFuncTestName<ChyokotovARunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(MatrixMultiplicate, ChyokotovARunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace chyokotov_a_dense_matrix_mul_foxs_algorithm
