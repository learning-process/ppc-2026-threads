#include "sannikov_i_integrals_rectangle_method/seq/include/ops_seq.hpp"

#include <cmath>
#include <cstddef>
#include <vector>

#include "sannikov_i_integrals_rectangle_method/common/include/common.hpp"

namespace sannikov_i_integrals_rectangle_method {

SannikovIIntegralsRectangleMethodSEQ::SannikovIIntegralsRectangleMethodSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool SannikovIIntegralsRectangleMethodSEQ::ValidationImpl() {
  const auto &[func, borders, n] = GetInput();
  if (borders.empty()) {
    return false;
  }
  for (const auto &[left_border, right_border] : borders) {
    if (!std::isfinite(left_border) || !std::isfinite(right_border)) {
      return false;
    }
    if (left_border >= right_border) {
      return false;
    }
  }

  return func && (n > 0) && (GetOutput() == 0.0);
}

bool SannikovIIntegralsRectangleMethodSEQ::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

namespace {
bool NextIndex(std::vector<int> &idx, std::size_t dim, int n) {
  for (std::size_t pos = 0; pos < dim; ++pos) {
    ++idx[pos];
    if (idx[pos] < n) {
      return true;
    }
    idx[pos] = 0;
  }
  return false;
}

}  // namespace

bool SannikovIIntegralsRectangleMethodSEQ::RunImpl() {
  const auto &[func, borders, n] = GetInput();
  const std::size_t dim = borders.size();

  std::vector<double> h(dim);
  double cell_v = 1.0;

  for (std::size_t i = 0; i < dim; ++i) {
    const double left_border = borders[i].first;
    const double right_border = borders[i].second;

    h[i] = (right_border - left_border) / static_cast<double>(n);
    if (!(h[i] > 0.0) || !std::isfinite(h[i])) {
      return false;
    }

    cell_v *= h[i];
  }

  std::vector<int> idx(dim, 0);
  std::vector<double> x(dim);

  double sum = 0.0;

  while (true) {
    for (std::size_t i = 0; i < dim; ++i) {
      const double left_border = borders[i].first;
      x[i] = left_border + ((static_cast<double>(idx[i]) + 0.5) * h[i]);
    }

    const double fx = func(x);
    if (!std::isfinite(fx)) {
      return false;
    }

    sum += fx;

    if (!NextIndex(idx, dim, n)) {
      break;
    }
  }

  GetOutput() = sum * cell_v;

  return std::isfinite(GetOutput());
}

bool SannikovIIntegralsRectangleMethodSEQ::PostProcessingImpl() {
  return std::isfinite(GetOutput());
}
}  // namespace sannikov_i_integrals_rectangle_method
