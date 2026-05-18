#include "kutergin_a_multidim_trapezoid/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <utility>
#include <vector>

#include "kutergin_a_multidim_trapezoid/common/include/common.hpp"

namespace kutergin_a_multidim_trapezoid {
namespace {

bool ValidateBorders(const std::vector<std::pair<double, double>> &borders) {
  return std::ranges::all_of(
      borders, [](const auto &p) { return std::isfinite(p.first) && std::isfinite(p.second) && (p.first < p.second); });
}

}  // namespace

KuterginAMultidimTrapezoidTBB::KuterginAMultidimTrapezoidTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool KuterginAMultidimTrapezoidTBB::ValidationImpl() {
  const auto &[func, borders, n] = GetInput();
  return func && n > 0 && !borders.empty() && ValidateBorders(borders);
}

bool KuterginAMultidimTrapezoidTBB::PreProcessingImpl() {
  local_input_ = GetInput();
  res_ = 0.0;
  return true;
}

bool KuterginAMultidimTrapezoidTBB::RunImpl() {
  const auto &[func, borders, n] = local_input_;
  const int dim = static_cast<int>(borders.size());

  std::vector<double> h(dim);
  double cell_volume = 1.0;
  for (int i = 0; i < dim; ++i) {
    h[i] = (borders[i].second - borders[i].first) / n;
    cell_volume *= h[i];
  }

  size_t total_points = 1;
  for (int i = 0; i < dim; ++i) {
    total_points *= (static_cast<size_t>(n) + 1);
  }

  res_ = tbb::parallel_reduce(tbb::blocked_range<size_t>(0, total_points), 0.0,
                              [=](const tbb::blocked_range<size_t> &r, double local_sum) {
    std::vector<double> point(dim);

    for (size_t i = r.begin(); i < r.end(); ++i) {
      double weight = 1.0;
      size_t temp_idx = i;

      for (int dim_idx = 0; dim_idx < dim; ++dim_idx) {
        int coord_idx = static_cast<int>(temp_idx % (n + 1));
        temp_idx /= (n + 1);

        point[dim_idx] = borders[dim_idx].first + (coord_idx * h[dim_idx]);

        if (coord_idx == 0 || coord_idx == n) {
          weight *= 0.5;
        }
      }
      local_sum += weight * func(point);
    }
    return local_sum;
  }, std::plus<>()) *
         cell_volume;

  return std::isfinite(res_);
}

bool KuterginAMultidimTrapezoidTBB::PostProcessingImpl() {
  GetOutput() = res_;
  return true;
}

}  // namespace kutergin_a_multidim_trapezoid
