#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <random>
#include <string>

#include "akhmetov_daniil_strassen_dense_double_seq/common/include/common.hpp"
#include "akhmetov_daniil_strassen_dense_double_seq/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"

namespace akhmetov_daniil_strassen_dense_double_seq {

namespace {

class AkhmetovDaniilRunFuncTests : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(test_param);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    InType input = GetTestInputData();

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

    constexpr double kEpsilon = 1e-7;
    for (size_t i = 0; i < n * n; ++i) {
      if (std::abs(output_data.at(i) - expected.at(i)) > kEpsilon) {
        return false;
      }
    }
    return true;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 protected:
  void SetUp() override {
    ppc::util::BaseRunFuncTests<InType, OutType, TestType>::SetUp();

    TestType n = std::get<2>(GetParam());
    input_data_.resize(1 + (2 * n * n));
    input_data_.at(0) = static_cast<double>(n);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(-10.0, 10.0);

    for (size_t i = 1; i < input_data_.size(); ++i) {
      input_data_.at(i) = dist(gen);
    }
  }

 private:
  InType input_data_;
};

TEST_P(AkhmetovDaniilRunFuncTests, StrassenTestFunctional) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 3> kTestParam = {64, 128, 256};

const auto kTestTasksList = ppc::util::AddFuncTask<AkhmetovDStrassenDenseDoubleSEQ, InType>(
    kTestParam, "tasks/akhmetov_daniil_strassen_dense_double_seq/settings.json");

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);
const auto kTestName = AkhmetovDaniilRunFuncTests::PrintFuncTestName<AkhmetovDaniilRunFuncTests>;

INSTANTIATE_TEST_SUITE_P(RunStrassenFuncTests, AkhmetovDaniilRunFuncTests, kGtestValues, kTestName);

}  // namespace
}  // namespace akhmetov_daniil_strassen_dense_double_seq
