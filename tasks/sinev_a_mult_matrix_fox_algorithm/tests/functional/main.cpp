#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "sinev_a_mult_matrix_fox_algorithm/common/include/common.hpp"
#include "sinev_a_mult_matrix_fox_algorithm/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace sinev_a_mult_matrix_fox_algorithm_seq {

class SinevARunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
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

    const double start = 1.0;
    const double step = 1.0;

    std::vector<double> a(size);
    std::vector<double> b(size);

    for (size_t i = 0; i < size; ++i) {
      a[i] = start + (static_cast<double>(i) * step);
    }

    for (size_t i = 0; i < size; ++i) {
      b[i] = start + (static_cast<double>(size - 1 - i) * step);
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

    const double epsilon = 1e-10;
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
      for (size_t j = 0; j < n; ++j) {
        double sum = 0.0;
        for (size_t k = 0; k < n; ++k) {
          sum += a[(i * n) + k] * b[(k * n) + j];
        }
        c[(i * n) + j] = sum;
      }
    }
  }

 private:
  InType input_data_;
  OutType expected_output_;
};

namespace {

TEST_P(SinevARunFuncTestsThreads, MatMulFoxAlg) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 7> kTestParams = {std::make_tuple(1, "size_1x1"), std::make_tuple(2, "size_2x2"),
                                             std::make_tuple(3, "size_3x3"), std::make_tuple(4, "size_4x4"),
                                             std::make_tuple(5, "size_5x5"), std::make_tuple(6, "size_6x6"),
                                             std::make_tuple(7, "size_7x7")};

const auto kTestTasksList = ppc::util::AddFuncTask<SinevAMultMatrixFoxAlgorithmSEQ, InType>(
    kTestParams, PPC_SETTINGS_sinev_a_mult_matrix_fox_algorithm);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = SinevARunFuncTestsThreads::PrintFuncTestName<SinevARunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(MatMulFoxAlg_ArithmeticProgression, SinevARunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace sinev_a_mult_matrix_fox_algorithm_seq
