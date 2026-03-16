#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <tuple>
#include <vector>

#include "cheremkhin_a_matr_mult_cannon_alg/common/include/common.hpp"
#include "cheremkhin_a_matr_mult_cannon_alg/seq/include/ops_seq.hpp"
#include "util/include/perf_test_util.hpp"

namespace cheremkhin_a_matr_mult_cannon_alg {

class CheremkhinAPerformanceTest : public ppc::util::BaseRunPerfTests<InType, OutType> {
 protected:
  void SetUp() override {
    std::size_t sz = 640;
    std::size_t size = sz * sz;

    std::vector<double> m1(sz * sz);
    std::vector<double> m2(sz * sz);

    for (std::size_t i = 0; i < size; i++) {
      m1[i] = static_cast<double>(i % 100);
      m2[i] = static_cast<double>((i + 1) % 100);
    }

    input_data_ = std::make_tuple(sz, m1, m2);

    std::vector<double> res(size, 0.0);
    MatMul(m1, m2, res, static_cast<int>(sz));
    expected_output_ = res;
  }

  bool CheckTestOutputData(OutType &output_data) final {
    const auto &expected = expected_output_;
    const auto &actual = output_data;

    if (expected.size() != actual.size()) {
      return false;
    }

    const double eps = 1e-7;
    for (std::size_t i = 0; i < expected.size(); ++i) {
      if (std::abs(expected[i] - actual[i]) > eps) {
        return false;
      }
    }
    return true;
  }

  static void MatMul(const std::vector<double> &m1, const std::vector<double> &m2, std::vector<double> &res, int n) {
    for (int i = 0; i < n; ++i) {
      for (int k = 0; k < n; ++k) {
        double tmp = m1[(i * n) + k];
        for (int j = 0; j < n; ++j) {
          res[(i * n) + j] += tmp * m2[(k * n) + j];
        }
      }
    }
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
  std::vector<double> expected_output_;
};

TEST_P(CheremkhinAPerformanceTest, RunPerformanceTest) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks = ppc::util::MakeAllPerfTasks<InType, CheremkhinAMatrMultCannonAlgSEQ>(
    PPC_SETTINGS_cheremkhin_a_matr_mult_cannon_alg);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = CheremkhinAPerformanceTest::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(SequentialPerformanceTests, CheremkhinAPerformanceTest, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace cheremkhin_a_matr_mult_cannon_alg
