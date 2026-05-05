#include "kutergin_a_multidim_trapezoid/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <ranges>
#include <tuple>
#include <utility>
#include <vector>

#include "kutergin_a_multidim_trapezoid/common/include/common.hpp"

namespace kutergin_a_multidim_trapezoid {

KuterginAMultidimTrapezoidOMP::KuterginAMultidimTrapezoidOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool KuterginAMultidimTrapezoidOMP::ValidationImpl() {
  const auto input = GetInput();
  auto func = std::get<0>(input);
  auto borders = std::get<1>(input);
  auto n = std::get<2>(input);

  if (!func || borders.empty() || n <= 0) {
    return false;
  }

  return std::all_of(borders.cbegin(), borders.cend(), [](const std::pair<double, double> &p) {
    return std::isfinite(p.first) && std::isfinite(p.second) && p.first < p.second;
  });
}

bool KuterginAMultidimTrapezoidOMP::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

bool KuterginAMultidimTrapezoidOMP::RunImpl() {
  const auto input = GetInput();
  auto func = std::get<0>(input);
  auto borders = std::get<1>(input);
  auto n = std::get<2>(input);

  const int dim = static_cast<int>(borders.size());

  std::vector<double> h(dim);
  for (int i = 0; i < dim; ++i) {
    h[i] = (borders[i].second - borders[i].first) / n;
  }

  int64_t total_points = 1;
  for (int i = 0; i < dim; ++i) {
    total_points *= (n + 1);
  }

  int prev_threads = omp_get_max_threads();
  int prev_dynamic = omp_get_dynamic();

  omp_set_dynamic(0);
  omp_set_num_threads(prev_threads);

  double global_sum = 0.0;

#pragma omp parallel default(none) shared(total_points, h, borders, func, n, dim) reduction(+ : global_sum)
  {
    std::vector<int> indices(dim);
    std::vector<double> point(dim);

#pragma omp for schedule(static)
    for (int64_t linear_idx = 0; linear_idx < total_points; ++linear_idx) {
      int64_t tmp = linear_idx;

      for (int axis = 0; axis < dim; ++axis) {
        indices[axis] = static_cast<int>(tmp % (n + 1));
        tmp /= (n + 1);
        point[axis] = borders[axis].first + (indices[axis] * h[axis]);
      }

      double weight = 1.0;
      for (int axis = 0; axis < dim; ++axis) {
        if (indices[axis] == 0 || indices[axis] == n) {
          weight *= 0.5;
        }
      }

      double val = func(point);
      if (std::isfinite(val)) {
        global_sum += weight * val;
      }
    }
  }

  omp_set_num_threads(prev_threads);
  omp_set_dynamic(prev_dynamic);

  double volume = 1.0;
  for (int i = 0; i < dim; ++i) {
    volume *= h[i];
  }

  GetOutput() = global_sum * volume;

  return std::isfinite(GetOutput());
}

bool KuterginAMultidimTrapezoidOMP::PostProcessingImpl() {
  return std::isfinite(GetOutput());
}

}  // namespace kutergin_a_multidim_trapezoid
