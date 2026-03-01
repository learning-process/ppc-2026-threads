#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <tuple>
#include <vector>

#include "sinev_a_mult_matrix_fox_algorithm/common/include/common.hpp"
#include "sinev_a_mult_matrix_fox_algorithm/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace sinev_a_mult_matrix_fox_algorithm_seq {

// ИСПРАВЛЕНО: Уникальное имя класса
class SinevAPerformanceTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;
  std::vector<double> expected_output_;

  void SetUp() override {
    size_t n = 400;

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
    const auto &expected = expected_output_;
    const auto &actual = output_data;

    if (expected.size() != actual.size()) {
      return false;
    }

    const double epsilon = 1e-7;
    for (size_t i = 0; i < expected.size(); ++i) {
      if (std::abs(expected[i] - actual[i]) > epsilon) {
        return false;
      }
    }
    return true;
  }

  static void ReferenceMultiply(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c,
                                size_t n) {
    for (size_t i = 0; i < n; ++i) {
      for (size_t k = 0; k < n; ++k) {
        double tmp = a[(i * n) + k];
        for (size_t j = 0; j < n; ++j) {
          c[(i * n) + j] += tmp * b[(k * n) + j];
        }
      }
    }
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(SinevAPerformanceTest, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, SinevAMultMatrixFoxAlgorithmSEQ>(
    PPC_SETTINGS_sinev_a_mult_matrix_fox_algorithm);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = SinevAPerformanceTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, SinevAPerformanceTest, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace sinev_a_mult_matrix_fox_algorithm_seq
