#include "sabirov_s_monte_carlo/tbb/include/ops_tbb.hpp"

#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/parallel_reduce.h"
#include "sabirov_s_monte_carlo/common/include/common.hpp"

namespace sabirov_s_monte_carlo {

SabirovSMonteCarloTBB::SabirovSMonteCarloTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool SabirovSMonteCarloTBB::ValidationImpl() {
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

bool SabirovSMonteCarloTBB::PreProcessingImpl() {
  const auto &in = GetInput();
  lower_ = in.lower;
  upper_ = in.upper;
  num_samples_ = in.num_samples;
  func_type_ = in.func_type;
  GetOutput() = 0.0;
  return true;
}

bool SabirovSMonteCarloTBB::RunImpl() {
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

  const double sum = tbb::parallel_reduce(tbb::blocked_range<int>(0, n_samples), 0.0,
                                          [&](const tbb::blocked_range<int> &r, double init) {
    double local = init;
    std::vector<double> point(static_cast<size_t>(dims));
    std::seed_seq seed{static_cast<uint32_t>(r.begin()), static_cast<uint32_t>(r.end()),
                       static_cast<uint32_t>(n_samples)};
    std::mt19937 gen(seed);
    for (int i = r.begin(); i != r.end(); ++i) {
      for (int j = 0; j < dims; ++j) {
        point[static_cast<size_t>(j)] = dists[static_cast<size_t>(j)](gen);
      }
      local += detail::EvaluateAt(ftype, point);
    }
    return local;
  }, [](double a, double b) { return a + b; });

  GetOutput() = volume * sum / static_cast<double>(num_samples_);
  return true;
}

bool SabirovSMonteCarloTBB::PostProcessingImpl() {
  return true;
}

}  // namespace sabirov_s_monte_carlo
