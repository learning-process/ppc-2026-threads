#include "dergynov_s_integrals_multistep_rectangle/omp/include/ops_omp.hpp"

#include <omp.h>

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

DergynovSIntegralsMultistepRectangleOMP::DergynovSIntegralsMultistepRectangleOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool DergynovSIntegralsMultistepRectangleOMP::ValidationImpl() {
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

bool DergynovSIntegralsMultistepRectangleOMP::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

bool DergynovSIntegralsMultistepRectangleOMP::RunImpl() {
  const auto &input = GetInput();
  const auto &func = std::get<0>(input);
  const auto &borders = std::get<1>(input);
  int n = std::get<2>(input);

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

  std::vector<double> local_sums(omp_get_max_threads(), 0.0);

  int error_flag = 0;

#pragma omp parallel default(none) shared(func, borders, h, dim, n, total_points, local_sums, error_flag, cell_volume)
  {
    int thread_id = omp_get_thread_num();
    double local_sum = 0.0;

#pragma omp for schedule(static)
    for (size_t linear_idx = 0; linear_idx < total_points; ++linear_idx) {
      if (error_flag != 0) {
        continue;
      }

      size_t tmp = linear_idx;
      std::vector<double> point(dim);

      for (int dimension = dim - 1; dimension >= 0; --dimension) {
        int idx_val = static_cast<int>(tmp % static_cast<size_t>(n));
        tmp /= static_cast<size_t>(n);

        point[dimension] = borders[dimension].first + ((static_cast<double>(idx_val) + 0.5) * h[dimension]);
      }

      double f_val = func(point);
      if (!std::isfinite(f_val)) {
#pragma omp atomic write
        error_flag = 1;
        continue;
      }
      local_sum += f_val;
    }

    local_sums[thread_id] = local_sum;
  }

  if (error_flag != 0) {
    return false;
  }

  double total_sum = 0.0;
  for (double s : local_sums) {
    total_sum += s;
  }

  GetOutput() = total_sum * cell_volume;
  return std::isfinite(GetOutput());
}

bool DergynovSIntegralsMultistepRectangleOMP::PostProcessingImpl() {
  return std::isfinite(GetOutput());
}

}  // namespace dergynov_s_integrals_multistep_rectangle
