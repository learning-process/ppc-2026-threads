#include "sabirov_s_monte_carlo/omp/include/ops_seq.hpp"

#include <cmath>
#include <cstddef>
#include <random>
#include <vector>

#include "sabirov_s_monte_carlo/common/include/common.hpp"

namespace sabirov_s_monte_carlo {

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

SabirovSMonteCarloOMP::SabirovSMonteCarloOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool SabirovSMonteCarloOMP::ValidationImpl() {
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

bool SabirovSMonteCarloOMP::PreProcessingImpl() {
  const auto &in = GetInput();
  lower_ = in.lower;
  upper_ = in.upper;
  num_samples_ = in.num_samples;
  func_type_ = in.func_type;
  GetOutput() = 0.0;
  return true;
}

bool SabirovSMonteCarloOMP::RunImpl() {
  auto dims = static_cast<int>(lower_.size());

  double volume = 1.0;
  for (int i = 0; i < dims; ++i) {
    volume *= (upper_[i] - lower_[i]);
  }

  double sum = 0.0;
  const std::vector<double> *const plower = &lower_;
  const std::vector<double> *const pupper = &upper_;
  const int n_samples = num_samples_;
  const FuncType ftype = func_type_;
#pragma omp parallel default(none) shared(plower, pupper, n_samples, ftype, dims, volume, sum)
  {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::vector<std::uniform_real_distribution<double>> dists;
    dists.reserve(static_cast<size_t>(dims));
    for (int i = 0; i < dims; ++i) {
      dists.emplace_back((*plower)[i], (*pupper)[i]);
    }
    std::vector<double> point(static_cast<size_t>(dims));

#pragma omp for reduction(+ : sum) schedule(static)
    for (int i = 0; i < n_samples; ++i) {
      for (int j = 0; j < dims; ++j) {
        point[j] = dists[j](gen);
      }
      sum += EvaluateAt(ftype, point);
    }
  }

  GetOutput() = volume * sum / static_cast<double>(num_samples_);
  return true;
}

bool SabirovSMonteCarloOMP::PostProcessingImpl() {
  return true;
}

}  // namespace sabirov_s_monte_carlo
