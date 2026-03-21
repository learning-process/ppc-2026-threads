#include "kutergin_a_multidim_trapezoid/omp/include/ops_omp.hpp"

#include <omp.h>

#include <cmath>
#include <vector>

#include "kutergin_a_multidim_trapezoid/common/include/common.hpp"

namespace kutergin_a_multidim_trapezoid {

KuterginAMultidimTrapezoidOMP::KuterginAMultidimTrapezoidOMP(const InType& in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool KuterginAMultidimTrapezoidOMP::ValidationImpl() {
  const auto& [func, borders, n] = GetInput();

  if (!func) return false;
  if (borders.empty()) return false;
  if (n <= 0) return false;

  for (const auto& [l, r] : borders) {
    if (!std::isfinite(l) || !std::isfinite(r) || l >= r) {
      return false;
    }
  }

  return true;
}

bool KuterginAMultidimTrapezoidOMP::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

bool KuterginAMultidimTrapezoidOMP::RunImpl() {
  const auto& [func, borders, n] = GetInput();
  const int dim = static_cast<int>(borders.size());

  // шаги
  std::vector<double> h(dim);
  for (int i = 0; i < dim; ++i) {
    h[i] = (borders[i].second - borders[i].first) / n;
  }

  // общее количество точек: (n+1)^dim
  long long total_points = 1;
  for (int i = 0; i < dim; ++i) {
    total_points *= (n + 1);
  }

  double global_sum = 0.0;

#pragma omp parallel reduction(+ : global_sum)
  {
    std::vector<int> idx(dim);
    std::vector<double> point(dim);

#pragma omp for schedule(static)
    for (long long linear = 0; linear < total_points; ++linear) {
      // раскодируем индекс
      long long tmp = linear;
      for (int d = 0; d < dim; ++d) {
        idx[d] = tmp % (n + 1);
        tmp /= (n + 1);
      }

      double weight = 1.0;

      for (int d = 0; d < dim; ++d) {
        point[d] = borders[d].first + idx[d] * h[d];

        if (idx[d] == 0 || idx[d] == n) {
          weight *= 0.5;
        }
      }

      double val = func(point);
      if (std::isfinite(val)) {
        global_sum += weight * val;
      }
    }
  }

  // объём ячейки
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