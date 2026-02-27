#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

// #include "korolev_k_matrix_mult/all/include/ops_all.hpp"
#include "korolev_k_matrix_mult/common/include/common.hpp"
// #include "korolev_k_matrix_mult/omp/include/ops_omp.hpp"
#include "korolev_k_matrix_mult/seq/include/ops_seq.hpp"
// #include "korolev_k_matrix_mult/stl/include/ops_stl.hpp"
// #include "korolev_k_matrix_mult/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"

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

MatrixInput CreateTestInput(size_t n) {
  MatrixInput in;
  in.n = n;
  in.A.resize(n * n);
  in.B.resize(n * n);
  for (size_t i = 0; i < n * n; ++i) {
    in.A[i] = static_cast<double>((i * 7 + 3) % 11 - 5);
    in.B[i] = static_cast<double>((i * 13 + 2) % 7 - 3);
  }
  return in;
}

}  // namespace

class KorolevKMatrixMultRunFuncTestsThreads : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(std::get<0>(test_param)) + "_" + std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    size_t n = std::get<0>(params);
    input_data_ = CreateTestInput(n);
    expected_ = NaiveMultiply(input_data_.A, input_data_.B, n);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    const size_t n = input_data_.n;
    constexpr double kTol = 1e-9;
    for (size_t i = 0; i < n * n; ++i) {
      if (std::fabs(output_data[i] - expected_[i]) > kTol) {
        return false;
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_{};
  std::vector<double> expected_;
};

namespace {

TEST_P(KorolevKMatrixMultRunFuncTestsThreads, MatrixMultiplyStrassen) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {std::make_tuple(4, "4"), std::make_tuple(8, "8"),
                                            std::make_tuple(16, "16")};

const auto kTestTasksList = std::tuple_cat(
    ppc::util::AddFuncTask<KorolevKMatrixMultALL, InType>(kTestParam, PPC_SETTINGS_korolev_k_matrix_mult),
    ppc::util::AddFuncTask<KorolevKMatrixMultOMP, InType>(kTestParam, PPC_SETTINGS_korolev_k_matrix_mult),
    ppc::util::AddFuncTask<KorolevKMatrixMultSEQ, InType>(kTestParam, PPC_SETTINGS_korolev_k_matrix_mult),
    ppc::util::AddFuncTask<KorolevKMatrixMultSTL, InType>(kTestParam, PPC_SETTINGS_korolev_k_matrix_mult),
    ppc::util::AddFuncTask<KorolevKMatrixMultTBB, InType>(kTestParam, PPC_SETTINGS_korolev_k_matrix_mult));

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName =
    KorolevKMatrixMultRunFuncTestsThreads::PrintFuncTestName<KorolevKMatrixMultRunFuncTestsThreads>;

INSTANTIATE_TEST_SUITE_P(MatrixMultiplyTests, KorolevKMatrixMultRunFuncTestsThreads, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace korolev_k_matrix_mult
