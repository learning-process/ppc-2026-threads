#include "kutergin_a_multidim_trapezoid/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

#include "kutergin_a_multidim_trapezoid/common/include/common.hpp"

namespace kutergin_a_multidim_trapezoid {
namespace {

bool ValidateBorders(const std::vector<std::pair<double, double>> &borders) {
  return std::ranges::all_of(
      borders, [](const auto &p) { return std::isfinite(p.first) && std::isfinite(p.second) && (p.first < p.second); });
}

bool NextIndex(std::vector<int> &idx, int dim, int max_index) {
  for (int pos = 0; pos < dim; ++pos) {
    ++idx[pos];
    if (idx[pos] <= max_index) {
      return true;
    }
    idx[pos] = 0;
  }
  return false;
}

}  // namespace

KuterginAMultidimTrapezoidSEQ::KuterginAMultidimTrapezoidSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool KuterginAMultidimTrapezoidSEQ::ValidationImpl() {
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

bool KuterginAMultidimTrapezoidSEQ::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

bool KuterginAMultidimTrapezoidSEQ::RunImpl() {
  const auto &[func, borders, n] = GetInput();
  const int dim = static_cast<int>(borders.size());

  std::vector<double> h(dim);
  double cell_volume = 1.0;

  for (int i = 0; i < dim; ++i) {
    const double left = borders[i].first;
    const double right = borders[i].second;
    h[i] = (right - left) / n;
    cell_volume *= h[i];
  }

  const int max_index = n;

  std::vector<int> idx(dim, 0);
  std::vector<double> point(dim);

  double sum = 0.0;

  while (true) {
    double weight = 1.0;

    for (int i = 0; i < dim; ++i) {
      point[i] = (borders[i].first + (idx[i] * h[i]));

      if ((idx[i] == 0) || (idx[i] == n)) {
        weight *= 0.5;
      }
    }

    double f_val = func(point);
    if (!std::isfinite(f_val)) {
      return false;
    }

    sum += weight * f_val;

    if (!NextIndex(idx, dim, max_index)) {
      break;
    }
  }

  GetOutput() = sum * cell_volume;
  return std::isfinite(GetOutput());
}

bool KuterginAMultidimTrapezoidSEQ::PostProcessingImpl() {
  return std::isfinite(GetOutput());
}

}  // namespace kutergin_a_multidim_trapezoid
