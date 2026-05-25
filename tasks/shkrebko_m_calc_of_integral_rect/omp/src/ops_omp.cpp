#include "shkrebko_m_calc_of_integral_rect/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "shkrebko_m_calc_of_integral_rect/common/include/common.hpp"
#include "util/include/util.hpp"

namespace shkrebko_m_calc_of_integral_rect {

ShkrebkoMCalcOfIntegralRectOMP::ShkrebkoMCalcOfIntegralRectOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool ShkrebkoMCalcOfIntegralRectOMP::ValidationImpl() {
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

bool ShkrebkoMCalcOfIntegralRectOMP::PreProcessingImpl() {
  local_input_ = GetInput();
  res_ = 0.0;
  return true;
}

bool ShkrebkoMCalcOfIntegralRectOMP::RunImpl() {
  const std::size_t dims = local_input_.limits.size();
  const auto &limits = local_input_.limits;
  const auto &n_steps = local_input_.n_steps;
  const auto &func = local_input_.func;

  std::int64_t total_points = 1;
  std::vector<double> h(dims);
  double cell_volume = 1.0;

  for (std::size_t i = 0; i < dims; ++i) {
    const double left = limits[i].first;
    const double right = limits[i].second;
    const int steps = n_steps[i];

    total_points *= static_cast<std::int64_t>(steps);
    h[i] = (right - left) / static_cast<double>(steps);
    cell_volume *= h[i];
  }

  double total_sum = 0.0;

  omp_set_num_threads(ppc::util::GetNumThreads());

#pragma omp parallel default(none) shared(total_points, h, limits, n_steps, func, dims) reduction(+ : total_sum)
  {
    std::vector<double> point(dims);

#pragma omp for
    for (std::int64_t idx = 0; idx < total_points; ++idx) {
      std::int64_t tmp = idx;

      for (int dim = static_cast<int>(dims) - 1; dim >= 0; --dim) {
        const int coord_index = static_cast<int>(tmp % n_steps[dim]);
        tmp /= n_steps[dim];

        point[dim] = limits[dim].first + ((static_cast<double>(coord_index) + 0.5) * h[dim]);
      }

      total_sum += func(point);
    }
  }

  res_ = total_sum * cell_volume;
  return true;
}
bool ShkrebkoMCalcOfIntegralRectOMP::PostProcessingImpl() {
  GetOutput() = res_;
  return true;
}

}  // namespace shkrebko_m_calc_of_integral_rect
