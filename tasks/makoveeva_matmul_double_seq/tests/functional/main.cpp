#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <random>
#include <string>
#include <tuple>
#include <vector>

#include "makoveeva_matmul_double_seq/common/include/common.hpp"
#include "makoveeva_matmul_double_seq/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace makoveeva_matmul_double_seq {

class MakoveevaRunFuncTestsSeq : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    size_t n = std::get<0>(test_param);
    std::string desc = std::get<1>(test_param);
    return desc + "_" + std::to_string(n) + "x" + std::to_string(n);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());

    size_t n = std::get<0>(params);
    size_t size = n * n;

    std::vector<double> a(size);
    std::vector<double> b(size);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(-50.0, 50.0);

    for (size_t i = 0; i < size; ++i) {
      a[i] = dist(gen);
      b[i] = dist(gen);
    }

    input_data_ = std::make_tuple(n, a, b);

    std::vector<double> expected(size, 0.0);
    ReferenceMultiply(a, b, expected, n);
    expected_output_ = expected;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    if (expected_output_.size() != output_data.size()) {
      return false;
    }

    const double epsilon = 1e-9;
    for (size_t i = 0; i < expected_output_.size(); ++i) {
      if (std::abs(expected_output_[i] - output_data[i]) > epsilon) {
        return false;
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

  static void ReferenceMultiply(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c,
                                size_t n) {
    for (size_t i = 0; i < n; ++i) {
      for (size_t k = 0; k < n; ++k) {
        const double a_ik = a[(i * n) + k];
        for (size_t j = 0; j < n; ++j) {
          c[(i * n) + j] += a_ik * b[(k * n) + j];
        }
      }
    }
  }

 private:
  InType input_data_;
  OutType expected_output_;
};

namespace {

TEST_P(MakoveevaRunFuncTestsSeq, FoxAlgorithmMatrixMultiplication) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 8> kTestParams = {std::make_tuple(1, "dimension_1"), std::make_tuple(2, "dimension_2"),
                                             std::make_tuple(3, "dimension_3"), std::make_tuple(4, "dimension_4"),
                                             std::make_tuple(6, "dimension_6"), std::make_tuple(8, "dimension_8"),
                                             std::make_tuple(9, "dimension_9"), std::make_tuple(12, "dimension_12")};

const auto kTestTasksList =
    ppc::util::AddFuncTask<MatmulDoubleSeqTask, InType>(kTestParams, PPC_SETTINGS_makoveeva_matmul_double_seq);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kFuncTestName = MakoveevaRunFuncTestsSeq::PrintFuncTestName<MakoveevaRunFuncTestsSeq>;

INSTANTIATE_TEST_SUITE_P(MatrixMultiplicationFoxSeq, MakoveevaRunFuncTestsSeq, kGtestValues, kFuncTestName);

}  // namespace

}  // namespace makoveeva_matmul_double_seq
