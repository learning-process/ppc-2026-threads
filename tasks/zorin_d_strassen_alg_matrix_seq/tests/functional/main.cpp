#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "util/include/func_test_util.hpp"
#include "zorin_d_strassen_alg_matrix_seq/common/include/common.hpp"
#include "zorin_d_strassen_alg_matrix_seq/seq/include/ops_seq.hpp"

namespace zorin_d_strassen_alg_matrix_seq {

namespace {

void NaiveMulRef(const std::vector<double> &A, const std::vector<double> &B, std::vector<double> &C, std::size_t n) {
  C.assign(n * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t k = 0; k < n; ++k) {
      const double aik = A[i * n + k];
      for (std::size_t j = 0; j < n; ++j) {
        C[i * n + j] += aik * B[k * n + j];
      }
    }
  }
}

std::vector<double> MakeMatrixA(std::size_t n) {
  std::vector<double> A(n * n);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      A[i * n + j] = std::sin(double(i + 1)) + std::cos(double(j + 1)) + double(i) * 0.1;
    }
  }
  return A;
}

std::vector<double> MakeMatrixB(std::size_t n) {
  std::vector<double> B(n * n);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      B[i * n + j] = std::cos(double(i + j + 1)) + double(j) * 0.05;
    }
  }
  return B;
}

class ZorinDRunFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &p) {
    return std::to_string(std::get<0>(p)) + "_" + std::get<1>(p);
  }

 protected:
  void SetUp() override {
    const auto params = std::get<std::size_t(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    const std::size_t n = static_cast<std::size_t>(std::get<0>(params));

    input_.n = n;
    input_.A = MakeMatrixA(n);
    input_.B = MakeMatrixB(n);

    NaiveMulRef(input_.A, input_.B, expected_, n);
  }

  InType GetTestInputData() final {
    return input_;
  }

  bool CheckTestOutputData(OutType &output) final {
    const double eps = 1e-9;
    if (output.size() != expected_.size()) {
      return false;
    }
    for (std::size_t i = 0; i < output.size(); ++i) {
      if (std::abs(output[i] - expected_[i]) > eps) {
        return false;
      }
    }
    return true;
  }

 private:
  InType input_{};
  std::vector<double> expected_;
};

TEST_P(ZorinDRunFuncTests, StrassenMatMul) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 7> kParams = {
    std::make_tuple(1, "n1"), std::make_tuple(2, "n2"), std::make_tuple(3, "n3"),   std::make_tuple(5, "n5"),
    std::make_tuple(8, "n8"), std::make_tuple(9, "n9"), std::make_tuple(16, "n16"),
};

const auto kTasks =
    ppc::util::AddFuncTask<ZorinDStrassenAlgMatrixSEQ, InType>(kParams, PPC_SETTINGS_zorin_d_strassen_alg_matrix_seq);

const auto kValues = ppc::util::ExpandToValues(kTasks);
const auto kName = ZorinDRunFuncTests::PrintFuncTestName<ZorinDRunFuncTests>;

INSTANTIATE_TEST_SUITE_P(StrassenMatrixTests, ZorinDRunFuncTests, kValues, kName);

}  // namespace

}  // namespace zorin_d_strassen_alg_matrix_seq
