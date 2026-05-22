#include <gtest/gtest.h>

#include <cmath>    // для std::abs
#include <cstddef>  // для size_t
#include <tuple>
#include <vector>

#include "makoveeva_matmul_double_seq/common/include/common.hpp"
#include "makoveeva_matmul_double_seq/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace makoveeva_matmul_double_seq {
namespace {

void ReferenceMultiply(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c, size_t n) {
  // Очищаем c перед вычислениями - используем assign вместо цикла
  c.assign(c.size(), 0.0);

  for (size_t i = 0; i < n; ++i) {
    for (size_t k = 0; k < n; ++k) {
      const double tmp = a[(i * n) + k];
      for (size_t j = 0; j < n; ++j) {
        c[(i * n) + j] += tmp * b[(k * n) + j];
      }
    }
  }
}

}  // namespace

class MatmulDoublePerformanceTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
  InType input_data_;
  std::vector<double> expected_output_;

 protected:
  void SetUp() override {
    const size_t n = 400;
    const size_t size = n * n;

    std::vector<double> a(size);
    std::vector<double> b(size);

    for (size_t i = 0; i < size; ++i) {
      a[i] = static_cast<double>(i + 1);
      b[i] = static_cast<double>(size - i);
    }

    input_data_ = std::make_tuple(n, a, b);

    expected_output_.resize(size);
    ReferenceMultiply(a, b, expected_output_, n);
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

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(MatmulDoublePerformanceTest, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, MatmulDoubleSeqTask>(PPC_SETTINGS_makoveeva_matmul_double_seq);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = MatmulDoublePerformanceTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, MatmulDoublePerformanceTest, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace makoveeva_matmul_double_seq
