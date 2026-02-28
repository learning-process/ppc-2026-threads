#include "dergynov_s_integrals_multistep_rectangle/seq/include/ops_seq.hpp"

#include <cmath>
#include <vector>

#include "dergynov_s_integrals_multistep_rectangle/common/include/common.hpp"

namespace dergynov_s_integrals_multistep_rectangle {
namespace {
bool ValidateBorders(const std::vector<std::pair<double, double>> &borders) {
  return std::ranges::all_of(
      borders, [](const auto &p) { return std::isfinite(p.first) && std::isfinite(p.second) && p.first < p.second; });
}

bool NextIndex(std::vector<int> &idx, int dim, int n) {
  for (int pos = 0; pos < dim; ++pos) {
    ++idx[pos];
    if (idx[pos] < n) {
      return true;
    }
    idx[pos] = 0;
  }
  return false;
}

}  // namespace

DergynovSIntegralsMultistepRectangleSEQ::DergynovSIntegralsMultistepRectangleSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool DergynovSIntegralsMultistepRectangleSEQ::ValidationImpl() {
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

bool DergynovSIntegralsMultistepRectangleSEQ::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

bool DergynovSIntegralsMultistepRectangleSEQ::RunImpl() {
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

  std::vector<int> idx(dim, 0);
  std::vector<double> point(dim);
  double sum = 0.0;

  while (true) {
    for (int i = 0; i < dim; ++i) {
      point[i] = borders[i].first + ((idx[i] + 0.5) * h[i]);
    }

    double f_val = func(point);
    if (!std::isfinite(f_val)) {
      return false;
    }
    sum += f_val;

    if (!NextIndex(idx, dim, n)) {
      break;
    }
  }

  GetOutput() = sum * cell_volume;
  return std::isfinite(GetOutput());
}

bool DergynovSIntegralsMultistepRectangleSEQ::PostProcessingImpl() {
  return std::isfinite(GetOutput());
}

}  // namespace dergynov_s_integrals_multistep_rectangle
