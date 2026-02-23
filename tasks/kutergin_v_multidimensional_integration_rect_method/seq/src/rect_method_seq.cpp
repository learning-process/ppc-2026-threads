#include "../include/rect_method_seq.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "../../common/include/common.hpp"

namespace kutergin_v_multidimensional_integration_rect_method {

RectMethodSequential::RectMethodSequential(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool RectMethodSequential::ValidationImpl() {
  const auto &input = GetInput();
  if (input.limits.size() != input.n_steps.size() || input.limits.empty()) {
    return false;
  }
  return std::ranges::all_of(input.n_steps, [](int n) { return n > 0; });
}

bool RectMethodSequential::PreProcessingImpl() {
  local_input_ = GetInput();
  res_ = 0.0;
  return true;
}

bool RectMethodSequential::RunImpl() {
  size_t dims = local_input_.limits.size();  // число размерностей пространства
  std::vector<double> coords(dims);          // создание вектора координат размером dims

  // вычисление гиперобъема
  double d_v = 1.0;
  for (size_t i = 0; i < dims; i++) {
    double h = (local_input_.limits[i].second - local_input_.limits[i].first) / local_input_.n_steps[i];
    d_v *= h;
  }

  res_ = CalculateRecursive(0, coords, d_v);
  return true;
}

bool RectMethodSequential::PostProcessingImpl() {
  GetOutput() = res_;
  return true;
}

// NOLINTNEXTLINE(misc-no-recursion)
double RectMethodSequential::CalculateRecursive(size_t dim, std::vector<double> coords, double d_v) {
  size_t dims = local_input_.limits.size();

  if (dim == dims) {
    return local_input_.func(coords) * d_v;
  }

  double total_sum = 0.0;
  double a = local_input_.limits[dim].first;
  double b = local_input_.limits[dim].second;
  int n = local_input_.n_steps[dim];
  double h = (b - a) / n;

  for (int i = 0; i < n; ++i) {
    coords[dim] = a + ((i + 0.5) * h);
    total_sum += CalculateRecursive(dim + 1, coords, d_v);
  }

  return total_sum;
}

}  // namespace kutergin_v_multidimensional_integration_rect_method
