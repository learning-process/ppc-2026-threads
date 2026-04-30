#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <random>
#include <string>

#include "akhmetov_daniil_strassen_dense_double_tbb/common/include/common.hpp"
#include "akhmetov_daniil_strassen_dense_double_tbb/tbb/include/ops_tbb.hpp"
#include "util/include/func_test_util.hpp"

namespace akhmetov_daniil_strassen_dense_double_tbb {

namespace {

class AkhmetovDaniilRunFuncTestsTBB : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::to_string(test_param);
  }

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

    const TestType n = std::get<2>(GetParam());
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

TEST_P(AkhmetovDaniilRunFuncTestsTBB, StrassenTestFunctional) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 10> kTestParam = {1, 2, 3, 5, 8, 9, 16, 64, 65, 257};

const auto kTestTasksList = ppc::util::AddFuncTask<AkhmetovDStrassenDenseDoubleTBB, InType>(
    kTestParam, "tasks/akhmetov_daniil_strassen_dense_double_tbb/settings.json");

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);
const auto kTestName = AkhmetovDaniilRunFuncTestsTBB::PrintFuncTestName<AkhmetovDaniilRunFuncTestsTBB>;

INSTANTIATE_TEST_SUITE_P(RunStrassenFuncTestsTBB, AkhmetovDaniilRunFuncTestsTBB, kGtestValues, kTestName);

TEST(AkhmetovDStrassenDenseDoubleTBBValidation, RejectsEmptyInput) {
  InType in;
  AkhmetovDStrassenDenseDoubleTBB task(in);
  EXPECT_FALSE(task.Validation());
}

TEST(AkhmetovDStrassenDenseDoubleTBBValidation, RejectsZeroSize) {
  InType in = {0.0};
  AkhmetovDStrassenDenseDoubleTBB task(in);
  EXPECT_FALSE(task.Validation());
}

TEST(AkhmetovDStrassenDenseDoubleTBBValidation, RejectsWrongInputSize) {
  constexpr int kN = 4;
  InType in(1 + (2 * kN * kN) - 1, 1.0);
  in.at(0) = static_cast<double>(kN);
  AkhmetovDStrassenDenseDoubleTBB task(in);
  EXPECT_FALSE(task.Validation());
}

}  // namespace
}  // namespace akhmetov_daniil_strassen_dense_double_tbb
