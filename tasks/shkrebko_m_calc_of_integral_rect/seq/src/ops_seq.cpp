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
  const auto &in = GetInput();

  if (!in.func) {
    return false;
  }
  if (in.limits.empty() || in.limits.size() != in.n_steps.size()) {
    return false;
  }
  if (!std::ranges::all_of(in.n_steps, [](int n) { return n > 0; })) {
    return false;
  }
  if (!std::ranges::all_of(in.limits, [](const auto &lim) { return lim.first < lim.second; })) {
    return false;
  }

  return true;
}

bool ShkrebkoMCalcOfIntegralRectSEQ::PreProcessingImpl() {
  return true;
}

bool ShkrebkoMCalcOfIntegralRectSEQ::RunImpl() {
  const auto &in = GetInput();
  const std::size_t dim = in.limits.size();

  std::vector<double> h(dim);
  double cell_volume = 1.0;
  for (std::size_t i = 0; i < dim; ++i) {
    const double left = in.limits[i].first;
    const double right = in.limits[i].second;
    const int steps = in.n_steps[i];
    h[i] = (right - left) / static_cast<double>(steps);
    cell_volume *= h[i];
  }

  std::vector<int> idx(dim, 0);
  std::vector<double> point(dim);
  double sum = 0.0;

  while (true) {
    for (std::size_t i = 0; i < dim; ++i) {
      point[i] = in.limits[i].first + ((static_cast<double>(idx[i]) + 0.5) * h[i]);
    }

    double f_val = in.func(point);
    if (!std::isfinite(f_val)) {
      return false;
    }
    sum += f_val;

    int level = static_cast<int>(dim) - 1;
    while (level >= 0) {
      if (++idx[level] < in.n_steps[level]) {
        break;
      }
      idx[level--] = 0;
    }
    if (level < 0) {
      break;
    }
  }

  GetOutput() = sum * cell_volume;
  return true;
}

bool ShkrebkoMCalcOfIntegralRectSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace shkrebko_m_calc_of_integral_rect
