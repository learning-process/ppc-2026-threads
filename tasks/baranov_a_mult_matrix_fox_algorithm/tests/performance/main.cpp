#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <tuple>
#include <vector>

#include "baranov_a_mult_matrix_fox_algorithm/common/include/common.hpp"
#include "baranov_a_mult_matrix_fox_algorithm/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace baranov_a_mult_matrix_fox_algorithm_seq {

class BaranovAPerfTest : public ppc::util::BaseRunPerfTests<baranov_a_mult_matrix_fox_algorithm::InType,
                                                            baranov_a_mult_matrix_fox_algorithm::OutType> {
  baranov_a_mult_matrix_fox_algorithm::InType input_data_;
  baranov_a_mult_matrix_fox_algorithm::OutType expected_output_;

  void SetUp() override {
    size_t n = 400;

    size_t size = n * n;

    std::vector<double> a(size);
    std::vector<double> b(size);

    for (size_t i = 0; i < size; ++i) {
      a[i] = std::pow(2.0, static_cast<double>(i % 20));
      b[i] = std::sqrt(static_cast<double>(i + 1));
    }

    input_data_ = std::make_tuple(n, a, b);

    std::vector<double> expected(size, 0.0);
    ReferenceMultiply(a, b, expected, n);
    expected_output_ = expected;
  }

  bool CheckTestOutputData(baranov_a_mult_matrix_fox_algorithm::OutType &output_data) final {
    const auto &expected = expected_output_;
    const auto &actual = output_data;

    if (expected.size() != actual.size()) {
      return false;
    }

    const double relative_epsilon = 1e-6;
    for (size_t i = 0; i < expected.size(); ++i) {
      double max_val = std::max(std::abs(expected[i]), std::abs(actual[i]));
      if (max_val > 1e-12) {
        double rel_error = std::abs(expected[i] - actual[i]) / max_val;
        if (rel_error > relative_epsilon) {
          return false;
        }
      } else {
        if (std::abs(expected[i] - actual[i]) > 1e-12) {
          return false;
        }
      }
    }
    return true;
  }

  static void ReferenceMultiply(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c,
                                size_t n) {
    for (size_t j = 0; j < n; ++j) {
      for (size_t i = 0; i < n; ++i) {
        double sum = 0.0;
        for (size_t k = 0; k < n; ++k) {
          sum += a[(i * n) + k] * b[(k * n) + j];
        }
        c[(i * n) + j] = sum;
      }
    }
  }

  baranov_a_mult_matrix_fox_algorithm::InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(BaranovAPerfTest, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<baranov_a_mult_matrix_fox_algorithm::InType, BaranovAMultMatrixFoxAlgorithmSEQ>(
        PPC_SETTINGS_baranov_a_mult_matrix_fox_algorithm);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = BaranovAPerfTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, BaranovAPerfTest, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace baranov_a_mult_matrix_fox_algorithm_seq
