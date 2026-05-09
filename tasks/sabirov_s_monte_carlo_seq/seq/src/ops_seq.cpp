#include "sabirov_s_monte_carlo_seq/seq/include/ops_seq.hpp"

#include <cmath>
#include <cstddef>
#include <random>
#include <vector>

#include "sabirov_s_monte_carlo_seq/common/include/common.hpp"

namespace sabirov_s_monte_carlo_seq {

namespace {

double EvalLinear(const std::vector<double> &point) {
  double s = 0.0;
  for (double x : point) {
    s += x;
  }
  return s;
}

double EvalSumCubes(const std::vector<double> &point) {
  double s = 0.0;
  for (double x : point) {
    s += x * x * x;
  }
  return s;
}

double EvalCosProduct(const std::vector<double> &point) {
  double p = 1.0;
  for (double x : point) {
    p *= std::cos(x);
  }
  return p;
}

double EvalExpNeg(const std::vector<double> &point) {
  double s = 0.0;
  for (double x : point) {
    s += x;
  }
  return std::exp(-s);
}

double EvalMixedPoly(const std::vector<double> &point) {
  double s = 0.0;
  for (double x : point) {
    s += (x * x) + x;
  }
  return s;
}

double EvalSinSum(const std::vector<double> &point) {
  double s = 0.0;
  for (double x : point) {
    s += std::sin(x);
  }
  return s;
}

double EvalSqrtSum(const std::vector<double> &point) {
  double s = 0.0;
  for (double x : point) {
    s += std::sqrt(x);
  }
  return s;
}

double EvalQuarticSum(const std::vector<double> &point) {
  double s = 0.0;
  for (double x : point) {
    s += x * x * x * x;
  }
  return s;
}

double EvaluateAt(FuncType func_type, const std::vector<double> &point) {
  switch (func_type) {
    case FuncType::kLinear:
      return EvalLinear(point);
    case FuncType::kSumCubes:
      return EvalSumCubes(point);
    case FuncType::kCosProduct:
      return EvalCosProduct(point);
    case FuncType::kExpNeg:
      return EvalExpNeg(point);
    case FuncType::kMixedPoly:
      return EvalMixedPoly(point);
    case FuncType::kSinSum:
      return EvalSinSum(point);
    case FuncType::kSqrtSum:
      return EvalSqrtSum(point);
    case FuncType::kQuarticSum:
      return EvalQuarticSum(point);
    default:
      return 0.0;
  }
}

}  // namespace

SabirovSMonteCarloSEQ::SabirovSMonteCarloSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool SabirovSMonteCarloSEQ::ValidationImpl() {
  const auto &in = GetInput();
  if (in.lower.size() != in.upper.size() || in.lower.empty()) {
    return false;
  }
  if (in.num_samples <= 0) {
    return false;
  }
  for (size_t i = 0; i < in.lower.size(); ++i) {
    if (in.lower[i] >= in.upper[i]) {
      return false;
    }
  }
  if (in.func_type < FuncType::kLinear || in.func_type > FuncType::kQuarticSum) {
    return false;
  }
  constexpr size_t kMaxDimensions = 10;
  return in.lower.size() <= kMaxDimensions;
}

bool SabirovSMonteCarloSEQ::PreProcessingImpl() {
  const auto &in = GetInput();
  lower_ = in.lower;
  upper_ = in.upper;
  num_samples_ = in.num_samples;
  func_type_ = in.func_type;
  GetOutput() = 0.0;
  return true;
}

bool SabirovSMonteCarloSEQ::RunImpl() {
  auto dims = static_cast<int>(lower_.size());

  std::random_device rd;
  std::mt19937 gen(rd());
  std::vector<std::uniform_real_distribution<double>> dists;
  dists.reserve(dims);
  for (int i = 0; i < dims; ++i) {
    dists.emplace_back(lower_[i], upper_[i]);
  }

  double volume = 1.0;
  for (int i = 0; i < dims; ++i) {
    volume *= (upper_[i] - lower_[i]);
  }

  std::vector<double> point(dims);
  double sum = 0.0;
  for (int i = 0; i < num_samples_; ++i) {
    for (int j = 0; j < dims; ++j) {
      point[j] = dists[j](gen);
    }
    sum += EvaluateAt(func_type_, point);
  }

  GetOutput() = volume * sum / static_cast<double>(num_samples_);
  return true;
}

bool SabirovSMonteCarloSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace sabirov_s_monte_carlo_seq
