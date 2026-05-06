#include "kutergin_a_multidim_trapezoid/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
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

  return std::ranges::all_of(borders, [](const std::pair<double, double> &p) {
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

  double global_sum = 0.0;

#pragma omp parallel default(none) reduction(+ : global_sum) firstprivate(dim, n) shared(h, borders, func, total_points)
  {
    std::vector<double> point(dim);
    std::vector<int> idx(dim);

#pragma omp for schedule(static)
    for (int64_t linear_idx = 0; linear_idx < total_points; ++linear_idx) {
      int64_t tmp = linear_idx;

      for (int dim_idx = 0; dim_idx < dim; ++dim_idx) {
        idx[dim_idx] = static_cast<int>(tmp % (n + 1));
        tmp /= (n + 1);

        point[dim_idx] = borders[dim_idx].first + (idx[dim_idx] * h[dim_idx]);
      }

      double weight = 1.0;
      for (int dim_idx = 0; dim_idx < dim; ++dim_idx) {
        if (idx[dim_idx] == 0 || idx[dim_idx] == n) {
          weight *= 0.5;
        }
      }

      double val = func(point);

      if (std::isfinite(val)) {
        global_sum += weight * val;
      }
    }
  }

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
