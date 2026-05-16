#include "dergynov_s_integrals_multistep_rectangle/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "dergynov_s_integrals_multistep_rectangle/common/include/common.hpp"

namespace dergynov_s_integrals_multistep_rectangle {
namespace {

bool ValidateBorders(const std::vector<std::pair<double, double>> &borders) {
  return std::ranges::all_of(borders, [](const auto &p) {
    const auto &[left, right] = p;
    return std::isfinite(left) && std::isfinite(right) && left < right;
  });
}

}  // namespace

DergynovSIntegralsMultistepRectangleTBB::DergynovSIntegralsMultistepRectangleTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool DergynovSIntegralsMultistepRectangleTBB::ValidationImpl() {
  const auto &[func, borders, n] = GetInput();

  if (!func) {
    return false;
  }
  if (n <= 0) {
    return false;
  }
  if (borders.empty()) {
    return false;
  }

  return ValidateBorders(borders);
}

bool DergynovSIntegralsMultistepRectangleTBB::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

bool DergynovSIntegralsMultistepRectangleTBB::RunImpl() {
  const auto &input = GetInput();
  const auto &func = std::get<0>(input);
  const auto &borders = std::get<1>(input);
  const int n = std::get<2>(input);
  const int dim = static_cast<int>(borders.size());

  std::vector<double> h(dim);
  double cell_volume = 1.0;

  for (int i = 0; i < dim; ++i) {
    const double left = borders[i].first;
    const double right = borders[i].second;
    h[i] = (right - left) / n;
    cell_volume *= h[i];
  }

  size_t total_points = 1;
  for (int i = 0; i < dim; ++i) {
    total_points *= n;
  }

  double total_sum = tbb::parallel_reduce(tbb::blocked_range<size_t>(0, total_points), 0.0,
                                          [&](const tbb::blocked_range<size_t> &range, double local_sum) {
    for (size_t linear_idx = range.begin(); linear_idx != range.end(); ++linear_idx) {
      size_t tmp = linear_idx;
      std::vector<double> point(dim);

      for (int dimension = dim - 1; dimension >= 0; --dimension) {
        int idx_val = static_cast<int>(tmp % static_cast<size_t>(n));
        tmp /= static_cast<size_t>(n);
        point[dimension] = borders[dimension].first + ((static_cast<double>(idx_val) + 0.5) * h[dimension]);
      }

      double f_val = func(point);
      if (!std::isfinite(f_val)) {
        return local_sum;
      }
      local_sum += f_val;
    }
    return local_sum;
  }, [](double a, double b) { return a + b; });

  GetOutput() = total_sum * cell_volume;
  return std::isfinite(GetOutput());
}

bool DergynovSIntegralsMultistepRectangleTBB::PostProcessingImpl() {
  return std::isfinite(GetOutput());
}

}  // namespace dergynov_s_integrals_multistep_rectangle
