#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "remizov_k_dense_matrix_multiplication_cannon_algorithm/common/include/common.hpp"
#include "remizov_k_dense_matrix_multiplication_cannon_algorithm/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace remizov_k_dense_matrix_multiplication_cannon_algorithm {

class RemizovKDenseMatrixMultiplicationCannonAlgorithmFuncTests
    : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<0>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());

    input_data_ = std::make_tuple(std::get<1>(params), std::get<2>(params), std::get<3>(params));
    expected_output_ = std::get<4>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (expected_output_.size() != output_data.size()) {
      return false;
    }

    if (expected_output_.empty()) {
      return true;
    }

    if (expected_output_[0].size() != output_data[0].size()) {
      return false;
    }

    double tolerance = 1e-10;
    for (size_t i = 0; i < expected_output_.size(); ++i) {
      for (size_t j = 0; j < expected_output_[0].size(); ++j) {
        if (std::abs(expected_output_[i][j] - output_data[i][j]) > tolerance) {
          return false;
        }
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

TEST_P(RemizovKDenseMatrixMultiplicationCannonAlgorithmFuncTests, MultiplicationMatrixBlockSchemeCannon) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 8> kTestCases = {
    std::make_tuple("test_1", 1, std::vector<std::vector<double>>{{2.0, 0.0}, {0.0, 2.0}},
                    std::vector<std::vector<double>>{{1.0, 2.0}, {3.0, 4.0}},
                    std::vector<std::vector<double>>{{2.0, 4.0}, {6.0, 8.0}}),

    std::make_tuple("test_2", 1, std::vector<std::vector<double>>{{1.0, 2.0}, {3.0, 4.0}},
                    std::vector<std::vector<double>>{{5.0, 6.0}, {7.0, 8.0}},
                    std::vector<std::vector<double>>{{19.0, 22.0}, {43.0, 50.0}}),

    std::make_tuple("test_3", 2, std::vector<std::vector<double>>(4, std::vector<double>(4, 1.0)),
                    std::vector<std::vector<double>>(4, std::vector<double>(4, 2.0)),
                    std::vector<std::vector<double>>(4, std::vector<double>(4, 8.0))),

    std::make_tuple("test_4", 2, std::vector<std::vector<double>>(4, std::vector<double>(4, 0.5)),
                    std::vector<std::vector<double>>(4, std::vector<double>(4, 0.5)),
                    std::vector<std::vector<double>>(4, std::vector<double>(4, 1.0))),

    std::make_tuple("test_5", 3, std::vector<std::vector<double>>{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}},
                    std::vector<std::vector<double>>{{9, 8, 7}, {6, 5, 4}, {3, 2, 1}},
                    std::vector<std::vector<double>>{{30, 24, 18}, {84, 69, 54}, {138, 114, 90}}),

    std::make_tuple("test_6", 3, std::vector<std::vector<double>>(6, std::vector<double>(6, 1.0)),
                    std::vector<std::vector<double>>(6, std::vector<double>(6, 1.0)),
                    std::vector<std::vector<double>>(6, std::vector<double>(6, 6.0))),

    std::make_tuple("test_7", 4, std::vector<std::vector<double>>(8, std::vector<double>(8, 1.0)),
                    std::vector<std::vector<double>>(8, std::vector<double>(8, 1.0)),
                    std::vector<std::vector<double>>(8, std::vector<double>(8, 8.0))),

    std::make_tuple("test_8", 2, std::vector<std::vector<double>>(4, std::vector<double>(4, 1.0)),
                    std::vector<std::vector<double>>(4, std::vector<double>(4, 1.0)),
                    std::vector<std::vector<double>>(4, std::vector<double>(4, 4.0)))};

const auto kTestTasksList =
    std::tuple_cat(ppc::util::AddFuncTask<RemizovKDenseMatrixMultiplicationCannonAlgorithm, InType>(
        kTestCases, PPC_SETTINGS_remizov_k_dense_matrix_multiplication_cannon_algorithm));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kTestNameFunc = RemizovKDenseMatrixMultiplicationCannonAlgorithmFuncTests::PrintFuncTestName<
    RemizovKDenseMatrixMultiplicationCannonAlgorithmFuncTests>;

INSTANTIATE_TEST_SUITE_P(CannonTests, RemizovKDenseMatrixMultiplicationCannonAlgorithmFuncTests, kGtestValues,
                         kTestNameFunc);

}  // namespace

}  // namespace remizov_k_dense_matrix_multiplication_cannon_algorithm
