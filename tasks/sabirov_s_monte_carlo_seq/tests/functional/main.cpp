#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <ostream>
#include <string>
#include <tuple>
#include <vector>

#include "sabirov_s_monte_carlo_seq/common/include/common.hpp"
#include "sabirov_s_monte_carlo_seq/seq/include/ops_seq.hpp"
#include "util/include/func_test_util.hpp"
#include "util/include/util.hpp"

namespace sabirov_s_monte_carlo_seq {

namespace {

std::ostream &operator<<(std::ostream &os, const MCInput &in) {
  os << "MCInput{func=" << static_cast<int>(in.func_type) << ",n=" << in.num_samples << ",dims=" << in.lower.size()
     << "}";
  return os;
}

double ExactLinear(const std::vector<double> &lo, const std::vector<double> &hi, double vol) {
  double mean = 0.0;
  for (size_t i = 0; i < lo.size(); ++i) {
    mean += (lo[i] + hi[i]) / 2.0;
  }
  return vol * mean;
}

double ExactSumCubes(const std::vector<double> &lo, const std::vector<double> &hi, double vol) {
  double mean = 0.0;
  for (size_t i = 0; i < lo.size(); ++i) {
    double b = hi[i];
    double a = lo[i];
    mean += (b * b * b + b * b * a + b * a * a + a * a * a) / 4.0;
  }
  return vol * mean;
}

double ExactCosProduct(const std::vector<double> &lo, const std::vector<double> &hi) {
  double prod = 1.0;
  for (size_t i = 0; i < lo.size(); ++i) {
    prod *= (std::sin(hi[i]) - std::sin(lo[i]));
  }
  return prod;
}

double ExactExpNeg(const std::vector<double> &lo, const std::vector<double> &hi) {
  double prod = 1.0;
  for (size_t i = 0; i < lo.size(); ++i) {
    prod *= (std::exp(-lo[i]) - std::exp(-hi[i]));
  }
  return prod;
}

double ExactMixedPoly(const std::vector<double> &lo, const std::vector<double> &hi, double vol) {
  double mean = 0.0;
  for (size_t i = 0; i < lo.size(); ++i) {
    double b = hi[i];
    double a = lo[i];
    mean += ((b * b + b * a + a * a) / 3.0) + ((b + a) / 2.0);
  }
  return vol * mean;
}

double ExactSinSum(const std::vector<double> &lo, const std::vector<double> &hi, double vol) {
  double mean = 0.0;
  for (size_t i = 0; i < lo.size(); ++i) {
    mean += (std::cos(lo[i]) - std::cos(hi[i])) / (hi[i] - lo[i]);
  }
  return vol * mean;
}

double ExactSqrtSum(const std::vector<double> &lo, const std::vector<double> &hi, double vol) {
  double mean = 0.0;
  for (size_t i = 0; i < lo.size(); ++i) {
    double b = hi[i];
    double a = lo[i];
    mean += (2.0 / 3.0) * (std::pow(b, 1.5) - std::pow(a, 1.5)) / (b - a);
  }
  return vol * mean;
}

double ExactQuarticSum(const std::vector<double> &lo, const std::vector<double> &hi, double vol) {
  double mean = 0.0;
  for (size_t i = 0; i < lo.size(); ++i) {
    double b = hi[i];
    double a = lo[i];
    mean += (b * b * b * b * b - a * a * a * a * a) / (5.0 * (b - a));
  }
  return vol * mean;
}

double ExactValue(const InType &in) {
  const auto &lo = in.lower;
  const auto &hi = in.upper;

  double vol = 1.0;
  for (size_t i = 0; i < lo.size(); ++i) {
    vol *= (hi[i] - lo[i]);
  }

  switch (in.func_type) {
    case FuncType::kLinear:
      return ExactLinear(lo, hi, vol);
    case FuncType::kSumCubes:
      return ExactSumCubes(lo, hi, vol);
    case FuncType::kCosProduct:
      return ExactCosProduct(lo, hi);
    case FuncType::kExpNeg:
      return ExactExpNeg(lo, hi);
    case FuncType::kMixedPoly:
      return ExactMixedPoly(lo, hi, vol);
    case FuncType::kSinSum:
      return ExactSinSum(lo, hi, vol);
    case FuncType::kSqrtSum:
      return ExactSqrtSum(lo, hi, vol);
    case FuncType::kQuarticSum:
      return ExactQuarticSum(lo, hi, vol);
    default:
      return 0.0;
  }
}

}  // namespace

class SabirovSMonteCarloFuncTestsSeq : public ppc::util::BaseRunFuncTests<InType, OutType, TestType> {
 public:
  static std::string PrintTestParam(const TestType &test_param) {
    return std::get<1>(test_param);
  }

 protected:
  void SetUp() override {
    TestType params = std::get<static_cast<std::size_t>(ppc::util::GTestParamIndex::kTestParams)>(GetParam());
    input_data_ = std::get<0>(params);
  }

  bool CheckTestOutputData(OutType &output_data) final {
    double tol = 5.0 / std::sqrt(static_cast<double>(input_data_.num_samples));
    return std::abs(output_data - ExactValue(input_data_)) <= tol;
  }

  InType GetTestInputData() final {
    return input_data_;
  }

 private:
  InType input_data_{};
};

namespace {

const double kHalfPi = std::acos(-1.0) / 2.0;

TEST_P(SabirovSMonteCarloFuncTestsSeq, MonteCarloIntegration) {
  ExecuteTest(GetParam());
}

const std::array<TestType, 13> kTestParam = {{
    std::make_tuple(MCInput{.lower = {0.0}, .upper = {1.0}, .num_samples = 10000, .func_type = FuncType::kLinear},
                    "kLinear_1D"),
    std::make_tuple(
        MCInput{.lower = {0.0, 0.0}, .upper = {1.0, 1.0}, .num_samples = 50000, .func_type = FuncType::kLinear},
        "kLinear_2D"),
    std::make_tuple(MCInput{.lower = {0.0}, .upper = {1.0}, .num_samples = 10000, .func_type = FuncType::kSumCubes},
                    "kSumCubes_1D"),
    std::make_tuple(
        MCInput{.lower = {0.0, 0.0}, .upper = {1.0, 1.0}, .num_samples = 50000, .func_type = FuncType::kSumCubes},
        "kSumCubes_2D"),
    std::make_tuple(
        MCInput{.lower = {0.0}, .upper = {kHalfPi}, .num_samples = 10000, .func_type = FuncType::kCosProduct},
        "kCosProduct_1D"),
    std::make_tuple(
        MCInput{
            .lower = {0.0, 0.0}, .upper = {kHalfPi, kHalfPi}, .num_samples = 50000, .func_type = FuncType::kCosProduct},
        "kCosProduct_2D"),
    std::make_tuple(MCInput{.lower = {0.0}, .upper = {1.0}, .num_samples = 10000, .func_type = FuncType::kExpNeg},
                    "kExpNeg_1D"),
    std::make_tuple(
        MCInput{.lower = {0.0, 0.0}, .upper = {1.0, 1.0}, .num_samples = 50000, .func_type = FuncType::kExpNeg},
        "kExpNeg_2D"),
    std::make_tuple(MCInput{.lower = {0.0}, .upper = {1.0}, .num_samples = 10000, .func_type = FuncType::kMixedPoly},
                    "kMixedPoly_1D"),
    std::make_tuple(MCInput{.lower = {0.0, 0.0, 0.0},
                            .upper = {1.0, 1.0, 1.0},
                            .num_samples = 100000,
                            .func_type = FuncType::kMixedPoly},
                    "kMixedPoly_3D"),
    std::make_tuple(MCInput{.lower = {0.0}, .upper = {kHalfPi}, .num_samples = 10000, .func_type = FuncType::kSinSum},
                    "kSinSum_1D"),
    std::make_tuple(MCInput{.lower = {0.0}, .upper = {1.0}, .num_samples = 10000, .func_type = FuncType::kSqrtSum},
                    "kSqrtSum_1D"),
    std::make_tuple(
        MCInput{.lower = {0.0, 0.0}, .upper = {1.0, 1.0}, .num_samples = 50000, .func_type = FuncType::kQuarticSum},
        "kQuarticSum_2D"),
}};

const auto kTestTasksList =
    ppc::util::AddFuncTask<SabirovSMonteCarloSEQ, InType>(kTestParam, PPC_SETTINGS_sabirov_s_monte_carlo_seq);

const auto kGtestValues = ppc::util::ExpandToValues(kTestTasksList);

const auto kPerfTestName = SabirovSMonteCarloFuncTestsSeq::PrintFuncTestName<SabirovSMonteCarloFuncTestsSeq>;

INSTANTIATE_TEST_SUITE_P(MonteCarloTests, SabirovSMonteCarloFuncTestsSeq, kGtestValues, kPerfTestName);

}  // namespace

}  // namespace sabirov_s_monte_carlo_seq
