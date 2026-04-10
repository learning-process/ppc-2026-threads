#include "shkrebko_m_calc_of_integral_rect/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

#include "shkrebko_m_calc_of_integral_rect/common/include/common.hpp"
#include "util/include/util.hpp"

namespace shkrebko_m_calc_of_integral_rect {

ShkrebkoMCalcOfIntegralRectTBB::ShkrebkoMCalcOfIntegralRectTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool ShkrebkoMCalcOfIntegralRectTBB::ValidationImpl() {
  const auto &input = GetInput();

  if (!input.func) {
    return false;
  }
  if (input.limits.size() != input.n_steps.size() || input.limits.empty()) {
    return false;
  }
  if (!std::ranges::all_of(input.n_steps, [](int n) { return n > 0; })) {
    return false;
  }
  if (!std::ranges::all_of(input.limits, [](const auto &lim) { return lim.first < lim.second; })) {
    return false;
  }

  return true;
}

bool ShkrebkoMCalcOfIntegralRectTBB::PreProcessingImpl() {
  local_input_ = GetInput();
  res_ = 0.0;
  return true;
}

bool ShkrebkoMCalcOfIntegralRectTBB::RunImpl() {
  const std::size_t dim = local_input_.limits.size();

  std::vector<double> h(dim);
  double cell_volume = 1.0;
  std::size_t total_points = 1;
  for (std::size_t i = 0; i < dim; ++i) {
    const double left = local_input_.limits[i].first;
    const double right = local_input_.limits[i].second;
    const int steps = local_input_.n_steps[i];
    h[i] = (right - left) / static_cast<double>(steps);
    cell_volume *= h[i];
    total_points *= static_cast<std::size_t>(steps);
  }

  int num_threads = ppc::util::GetNumThreads();
  tbb::global_control global_limit(tbb::global_control::max_allowed_parallelism, num_threads);

  double total_sum = tbb::parallel_reduce(tbb::blocked_range<std::size_t>(0, total_points), 0.0,
                                          [&](const tbb::blocked_range<std::size_t> &range, double partial_sum) {
    return partial_sum + ComputeBlockSum(range.begin(), range.end(), h);
  }, std::plus<>());

  res_ = total_sum * cell_volume;
  return true;
}

bool ShkrebkoMCalcOfIntegralRectTBB::PostProcessingImpl() {
  GetOutput() = res_;
  return true;
}

double ShkrebkoMCalcOfIntegralRectTBB::ComputeBlockSum(std::size_t start_idx, std::size_t end_idx,
                                                       const std::vector<double> &h) const {
  if (start_idx >= end_idx) {
    return 0.0;
  }

  const std::size_t dim = local_input_.limits.size();
  const auto &limits = local_input_.limits;
  const auto &n_steps = local_input_.n_steps;

  std::vector<int> indices(dim);
  std::size_t temp = start_idx;
  for (int idx = static_cast<int>(dim) - 1; idx >= 0; --idx) {
    {
      indices[idx] = static_cast<int>(temp % static_cast<std::size_t>(n_steps[idx]));
      temp /= static_cast<std::size_t>(n_steps[idx]);
    }
  }
  double block_sum = 0.0;
  std::vector<double> point(dim);

  for (std::size_t iter = 0; iter < (end_idx - start_idx); ++iter) {
    for (std::size_t dim_idx = 0; dim_idx < dim; ++dim_idx) {
      point[dim_idx] = limits[dim_idx].first + ((static_cast<double>(indices[dim_idx]) + 0.5) * h[dim_idx]);
    }
    block_sum += local_input_.func(point);

    int level = static_cast<int>(dim) - 1;
    while (level >= 0) {
      indices[level]++;
      if (indices[level] < n_steps[level]) {
        break;
      }
      indices[level] = 0;
      level--;
    }
  }

  return block_sum;
}

}  // namespace shkrebko_m_calc_of_integral_rect
