#include "kiselev_i_trapezoidal_method_for_multidimensional_integrals/seq/include/ops_seq.hpp"

#include <cmath>
#include <vector>

#include "kiselev_i_trapezoidal_method_for_multidimensional_integrals/common/include/common.hpp"

namespace kiselev_i_trapezoidal_method_for_multidimensional_integrals {

KiselevITestTaskSEQ::KiselevITestTaskSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool KiselevITestTaskSEQ::ValidationImpl() {
  return true;
}

bool KiselevITestTaskSEQ::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

double KiselevITestTaskSEQ::FunctionTypeChoose(int type_x, double x, double y) {
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

double KiselevITestTaskSEQ::ComputeIntegral(const std::vector<int> &steps) {
  double result = 0.0;

  double hx = (GetInput().right_bounds[0] - GetInput().left_bounds[0]) / steps[0];
  double hy = (GetInput().right_bounds[1] - GetInput().left_bounds[1]) / steps[1];

  for (int i = 0; i <= steps[0]; i++) {
    double x = GetInput().left_bounds[0] + (i * hx);
    double wx = (i == 0 || i == steps[0]) ? 0.5 : 1.0;

    for (int j = 0; j <= steps[1]; j++) {
      double y = GetInput().left_bounds[1] + (j * hy);
      double wy = (j == 0 || j == steps[1]) ? 0.5 : 1.0;

      result += wx * wy * FunctionTypeChoose(GetInput().type_function, x, y);
    }
  }

  return result * hx * hy;
}

bool KiselevITestTaskSEQ::RunImpl() {
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

bool KiselevITestTaskSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace kiselev_i_trapezoidal_method_for_multidimensional_integrals
