#include "shkrebko_m_calc_of_integral_rect/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "shkrebko_m_calc_of_integral_rect/common/include/common.hpp"

namespace shkrebko_m_calc_of_integral_rect {

ShkrebkoMCalcOfIntegralRectSEQ::ShkrebkoMCalcOfIntegralRectSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool ShkrebkoMCalcOfIntegralRectSEQ::ValidationImpl() {
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

bool ShkrebkoMCalcOfIntegralRectSEQ::PreProcessingImpl() {
  local_input_ = GetInput();
  res_ = 0.0;
  return true;
}

bool ShkrebkoMCalcOfIntegralRectSEQ::RunImpl() {
  const std::size_t dim = local_input_.limits.size();
  std::vector<double> h(dim);
  double cell_volume = 1.0;
  for (std::size_t i = 0; i < dim; ++i) {
    const double left = local_input_.limits[i].first;
    const double right = local_input_.limits[i].second;
    const int steps = local_input_.n_steps[i];
    h[i] = (right - left) / static_cast<double>(steps);
    cell_volume *= h[i];
  }

  std::vector<int> idx(dim, 0);
  std::vector<double> point(dim);
  double sum = 0.0;

  while (true) {
    for (std::size_t i = 0; i < dim; ++i) {
      point[i] = local_input_.limits[i].first + ((static_cast<double>(idx[i]) + 0.5) * h[i]);
    }

    double f_val = local_input_.func(point);
    if (!std::isfinite(f_val)) {
      return false;
    }
    sum += f_val;

    int level = static_cast<int>(dim) - 1;
    while (level >= 0) {
      idx[level]++;
      if (idx[level] < local_input_.n_steps[level]) {
        break;
      }
      idx[level] = 0;
      level--;
    }
    if (level < 0) {
      break;
    }
  }

  res_ = sum * cell_volume;
  return true;
}

bool ShkrebkoMCalcOfIntegralRectSEQ::PostProcessingImpl() {
  GetOutput() = res_;
  return true;
}

}  // namespace shkrebko_m_calc_of_integral_rect
