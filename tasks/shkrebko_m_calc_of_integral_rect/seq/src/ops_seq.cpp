#include "shkrebko_m_calc_of_integral_rect/seq/include/ops_seq.hpp"

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

  if (input.limits.empty()) {
    return false;
  }

  if (input.steps <= 0) {
    return false;
  }

  for (std::size_t i = 0; i < input.limits.size(); ++i) {
    if (input.limits[i].first >= input.limits[i].second) {
      return false;
    }
  }
  return true;
}

bool ShkrebkoMCalcOfIntegralRectSEQ::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

bool ShkrebkoMCalcOfIntegralRectSEQ::RunImpl() {
  const auto &input = GetInput();
  const std::size_t dim = input.limits.size();
  const int n = input.steps;

  if (dim == 0) {
    GetOutput() = 0.0;
    return true;
  }

  std::vector<double> h(dim);
  double cell_volume = 1.0;
  for (std::size_t i = 0; i < dim; ++i) {
    h[i] = (input.limits[i].second - input.limits[i].first) / static_cast<double>(n);
    cell_volume *= h[i];
  }

  long long total_cells = 1;
  for (std::size_t i = 0; i < dim; ++i) {
    total_cells *= n;
  }

  double sum = 0.0;
  std::vector<int> idx(dim);

  for (long long cell = 0; cell < total_cells; ++cell) {
    long long tmp = cell;
    for (std::size_t i = 0; i < dim; ++i) {
      idx[i] = static_cast<int>(tmp % n);
      tmp /= n;
    }

    std::vector<double> point(dim);
    for (std::size_t i = 0; i < dim; ++i) {
      point[i] = input.limits[i].first + ((static_cast<double>(idx[i]) + 0.5) * h[i]);
    }

    double f_val = input.func(point);
    if (!std::isfinite(f_val)) {
      return false;
    }
    sum += f_val;
  }

  GetOutput() = sum * cell_volume;
  return true;
}

bool ShkrebkoMCalcOfIntegralRectSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace shkrebko_m_calc_of_integral_rect
