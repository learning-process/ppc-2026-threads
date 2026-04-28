#include "kiselev_i_trapezoidal_method_for_multidimensional_integrals/omp/include/ops_omp.hpp"

#include <omp.h>

#include <cmath>
#include <vector>

#include "kiselev_i_trapezoidal_method_for_multidimensional_integrals/common/include/common.hpp"

namespace kiselev_i_trapezoidal_method_for_multidimensional_integrals {

KiselevITestTaskOMP::KiselevITestTaskOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool KiselevITestTaskOMP::ValidationImpl() {
  return true;
}

bool KiselevITestTaskOMP::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

double KiselevITestTaskOMP::FunctionTypeChoose(int type_x, double x, double y) {
  switch (type_x) {
    case 0:
      return (x * x) + (y * y);
    case 1:
      return std::sin(x) * std::cos(y);
    case 2:
      return std::sin(x) + std::cos(y);
    case 3:
      return std::exp(x + y);
    default:
      return x + y;
  }
}

double KiselevITestTaskOMP::ComputeIntegral(const std::vector<int> &steps) {
  double result = 0.0;

  double hx = (GetInput().right_bounds[0] - GetInput().left_bounds[0]) / steps[0];
  double hy = (GetInput().right_bounds[1] - GetInput().left_bounds[1]) / steps[1];

  const double x0 = GetInput().left_bounds[0];
  const double y0 = GetInput().left_bounds[1];

  int nx = steps[0];
  int ny = steps[1];

  const int func_type = GetInput().type_function;

#pragma omp parallel for reduction(+ : result) default(none) shared(nx, ny, hx, hy, x0, y0, func_type)
  for (int i = 0; i <= nx; i++) {
    for (int j = 0; j <= ny; j++) {
      const double x = x0 + (i * hx);
      const double y = y0 + (j * hy);

      const double wx = (i == 0 || i == nx) ? 0.5 : 1.0;
      const double wy = (j == 0 || j == ny) ? 0.5 : 1.0;

      result += wx * wy * FunctionTypeChoose(func_type, x, y);
    }
  }

  return result * hx * hy;
}

bool KiselevITestTaskOMP::RunImpl() {
  std::vector<int> steps = GetInput().step_n_size;
  double epsilon = GetInput().epsilon;

  const auto &in = GetInput();
  if (in.left_bounds.size() != 2 || in.right_bounds.size() != 2 || in.step_n_size.size() != 2) {
    GetOutput() = 0.0;
    return true;
  }
  if (epsilon <= 0.0) {
    GetOutput() = ComputeIntegral(steps);
    return true;
  }

  double prev = ComputeIntegral(steps);
  double current = prev;

  int iter = 0;
  const int max_iter = 1;  // for time_limit

  while (iter < max_iter) {
    for (auto &s : steps) {
      s *= 2;
    }

    current = ComputeIntegral(steps);

    if (std::abs(current - prev) < epsilon) {
      break;
    }

    prev = current;
    iter++;
  }

  GetOutput() = current;
  return true;
}

bool KiselevITestTaskOMP::PostProcessingImpl() {
  return true;
}

}  // namespace kiselev_i_trapezoidal_method_for_multidimensional_integrals
