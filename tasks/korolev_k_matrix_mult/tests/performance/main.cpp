#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include "korolev_k_matrix_mult/all/include/ops_all.hpp"
#include "korolev_k_matrix_mult/common/include/common.hpp"
#include "korolev_k_matrix_mult/omp/include/ops_omp.hpp"
#include "korolev_k_matrix_mult/seq/include/ops_seq.hpp"
#include "korolev_k_matrix_mult/stl/include/ops_stl.hpp"
#include "korolev_k_matrix_mult/tbb/include/ops_tbb.hpp"
#include "util/include/perf_test_util.hpp"

namespace korolev_k_matrix_mult {

namespace {

std::vector<double> NaiveMultiply(const std::vector<double> &A, const std::vector<double> &B, size_t n) {
  std::vector<double> C(n * n, 0.0);
  for (size_t i = 0; i < n; ++i) {
    for (size_t k = 0; k < n; ++k) {
      double a_ik = A[i * n + k];
      for (size_t j = 0; j < n; ++j) {
        C[i * n + j] += a_ik * B[k * n + j];
      }
    }
  }
  return C;
}

}  // namespace

class KorolevKMatrixMultRunPerfTestThreads : public ppc::util::BaseRunPerfTests<InType, OutType> {
  static constexpr size_t kMatrixSize = 256;
  InType input_data_{};

  void SetUp() override {
    input_data_.n = kMatrixSize;
    input_data_.A.resize(kMatrixSize * kMatrixSize);
    input_data_.B.resize(kMatrixSize * kMatrixSize);
    for (size_t i = 0; i < kMatrixSize * kMatrixSize; ++i) {
      input_data_.A[i] = static_cast<double>((i * 7 + 3) % 11 - 5);
      input_data_.B[i] = static_cast<double>((i * 13 + 2) % 7 - 3);
    }
  }

  bool CheckTestOutputData(OutType &output_data) final {
    auto expected = NaiveMultiply(input_data_.A, input_data_.B, kMatrixSize);
    constexpr double kTol = 1e-9;
    for (size_t i = 0; i < kMatrixSize * kMatrixSize; ++i) {
      if (std::fabs(output_data[i] - expected[i]) > kTol) {
        return false;
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }
};

TEST_P(KorolevKMatrixMultRunPerfTestThreads, RunPerfModes) {
  ExecuteTest(GetParam());
}

namespace {

const auto kAllPerfTasks =
    ppc::util::MakeAllPerfTasks<InType, KorolevKMatrixMultALL, KorolevKMatrixMultOMP, KorolevKMatrixMultSEQ,
                                KorolevKMatrixMultSTL, KorolevKMatrixMultTBB>(PPC_SETTINGS_korolev_k_matrix_mult);

const auto kGtestValues = ppc::util::TupleToGTestValues(kAllPerfTasks);

const auto kPerfTestName = KorolevKMatrixMultRunPerfTestThreads::CustomPerfTestName;

INSTANTIATE_TEST_SUITE_P(RunModeTests, KorolevKMatrixMultRunPerfTestThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace korolev_k_matrix_mult
