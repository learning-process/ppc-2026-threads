#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <random>
#include <string>

#include "akhmetov_daniil_strassen_dense_double_stl/common/include/common.hpp"
#include "akhmetov_daniil_strassen_dense_double_stl/stl/include/ops_stl.hpp"
#include "util/include/func_test_util.hpp"

namespace akhmetov_daniil_strassen_dense_double_stl {

namespace {

class AkhmetovDaniilRunFuncTestsSTL : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(test_param);
  }

 protected:
  void SetUp() override {
    ppc::util::BaseRunFuncTests<InType, OutType, TestType>::SetUp();

    const TestType n = std::get<2>(GetParam());
    input_data_.resize(1 + (2 * n * n));
    input_data_.at(0) = static_cast<double>(n);

    // Deterministic seed to avoid flaky CI failures due to random data.
    std::mt19937 gen(static_cast<std::mt19937::result_type>(0xC0FFEEu + static_cast<unsigned>(n)));
    std::uniform_real_distribution<double> dist(-10.0, 10.0);

    for (size_t i = 1; i < input_data_.size(); ++i) {
      input_data_.at(i) = dist(gen);
    }
  }

 public:
  bool CheckTestOutputData(OutType &output_data) final {
    const InType input = GetTestInputData();

    const size_t n = format::GetN(input);
    const Matrix a = format::GetA(input);
    const Matrix b = format::GetB(input);

    Matrix expected(n * n, 0.0);
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j) {
        double sum = 0.0;
        for (size_t k = 0; k < n; ++k) {
          sum += a.at((i * n) + k) * b.at((k * n) + j);
        }
        expected.at((i * n) + j) = sum;
      }
    }

    constexpr double kEpsilon = 1e-6;
    for (size_t i = 0; i < n * n; ++i) {
      const double exp = expected.at(i);
      const double diff = std::abs(output_data.at(i) - exp);
      const double tol = kEpsilon * (1.0 + std::abs(exp));
      if (diff > tol) {
        return false;
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_;
};

TEST_P(AkhmetovDaniilRunFuncTestsSTL, StrassenTestFunctional) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 10> kTestParam = {1, 2, 3, 5, 8, 9, 16, 64, 65, 257};

const auto kTestTasksList = ppc::util::AddFuncTask<AkhmetovDStrassenDenseDoubleSTL, InType>(
    kTestParam, "tasks/akhmetov_daniil_strassen_dense_double_stl/settings.json");

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);
const auto kTestName = AkhmetovDaniilRunFuncTestsSTL::PrintFuncTestName<AkhmetovDaniilRunFuncTestsSTL>;

INSTANTIATE_TEST_SUITE_P(RunStrassenFuncTestsSTL, AkhmetovDaniilRunFuncTestsSTL, kGtestValues, kTestName);

}  // namespace
}  // namespace akhmetov_daniil_strassen_dense_double_stl
