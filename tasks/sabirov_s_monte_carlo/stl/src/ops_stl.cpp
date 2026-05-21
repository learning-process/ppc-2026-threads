#include "sabirov_s_monte_carlo/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>
#include <thread>
#include <vector>

#include "sabirov_s_monte_carlo/common/include/common.hpp"
#include "util/include/util.hpp"

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

SabirovSMonteCarloSTL::SabirovSMonteCarloSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool SabirovSMonteCarloSTL::ValidationImpl() {
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

bool SabirovSMonteCarloSTL::PreProcessingImpl() {
  const auto &in = GetInput();
  lower_ = in.lower;
  upper_ = in.upper;
  num_samples_ = in.num_samples;
  func_type_ = in.func_type;
  GetOutput() = 0.0;
  return true;
}

bool SabirovSMonteCarloSTL::RunImpl() {
  const int dims = static_cast<int>(lower_.size());

  std::vector<std::uniform_real_distribution<double>> dists;
  dists.reserve(static_cast<size_t>(dims));
  for (int j = 0; j < dims; ++j) {
    dists.emplace_back(lower_[j], upper_[j]);
  }

  double volume = 1.0;
  for (int j = 0; j < dims; ++j) {
    volume *= (upper_[j] - lower_[j]);
  }

  const FuncType ftype = func_type_;
  const int n_samples = num_samples_;
  const int num_threads = std::max(1, ppc::util::GetNumThreads());

  std::vector<double> partial_sums(static_cast<size_t>(num_threads), 0.0);
  std::vector<std::thread> threads;
  threads.reserve(static_cast<size_t>(num_threads));

  for (int ti = 0; ti < num_threads; ++ti) {
    const int start = (ti * n_samples) / num_threads;
    const int end = ((ti + 1) * n_samples) / num_threads;
    threads.emplace_back([&, ti, start, end, dims, ftype]() {
      if (start >= end) {
        return;
      }
      std::vector<double> point(static_cast<size_t>(dims));
      std::seed_seq seed{static_cast<uint32_t>(ti), static_cast<uint32_t>(start), static_cast<uint32_t>(end),
                         static_cast<uint32_t>(n_samples)};
      std::mt19937 gen(seed);
      double local = 0.0;
      for (int i = start; i < end; ++i) {
        for (int j = 0; j < dims; ++j) {
          point[static_cast<size_t>(j)] = dists[static_cast<size_t>(j)](gen);
        }
        local += EvaluateAt(ftype, point);
      }
      partial_sums[static_cast<size_t>(ti)] = local;
    });
  }

  for (auto &th : threads) {
    th.join();
  }

  double sum = 0.0;
  for (double part : partial_sums) {
    sum += part;
  }

  GetOutput() = volume * sum / static_cast<double>(num_samples_);
  return true;
}

bool SabirovSMonteCarloSTL::PostProcessingImpl() {
  return true;
}

}  // namespace sabirov_s_monte_carlo
