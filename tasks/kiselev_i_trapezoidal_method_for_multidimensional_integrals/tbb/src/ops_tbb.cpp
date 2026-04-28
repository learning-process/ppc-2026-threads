#include "kiselev_i_trapezoidal_method_for_multidimensional_integrals/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range2d.h>
#include <tbb/parallel_reduce.h>

#include <cmath>
#include <functional>
#include <vector>

#include "kiselev_i_trapezoidal_method_for_multidimensional_integrals/common/include/common.hpp"

namespace kiselev_i_trapezoidal_method_for_multidimensional_integrals {

KiselevITestTaskTBB::KiselevITestTaskTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool KiselevITestTaskTBB::ValidationImpl() {
  return true;
}

bool KiselevITestTaskTBB::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

double KiselevITestTaskTBB::FunctionTypeChoose(int type_x, double x, double y) {
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

double KiselevITestTaskTBB::ComputeIntegral(const std::vector<int> &steps) {
  const auto &in = GetInput();

  double hx = (in.right_bounds[0] - in.left_bounds[0]) / steps[0];
  double hy = (in.right_bounds[1] - in.left_bounds[1]) / steps[1];
  using Range2d = tbb::blocked_range2d<int, int>;

  double result = tbb::parallel_reduce(Range2d(0, steps[0] + 1, 32, 0, steps[1] + 1, 32), 0.0,
                                       [&](const Range2d &r, double local_sum) {
    for (int i = r.rows().begin(); i < r.rows().end(); ++i) {
      double x = in.left_bounds[0] + (i * hx);
      double wx = (i == 0 || i == steps[0]) ? 0.5 : 1.0;

      for (int j = r.cols().begin(); j < r.cols().end(); ++j) {
        double y = in.left_bounds[1] + (j * hy);
        double wy = (j == 0 || j == steps[1]) ? 0.5 : 1.0;

        local_sum += wx * wy * FunctionTypeChoose(in.type_function, x, y);
      }
    }
    return local_sum;
  }, std::plus<>());

  return result * hx * hy;
}

bool KiselevITestTaskTBB::RunImpl() {
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
  const int max_iter = 10;

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

bool KiselevITestTaskTBB::PostProcessingImpl() {
  return true;
}

}  // namespace kiselev_i_trapezoidal_method_for_multidimensional_integrals
