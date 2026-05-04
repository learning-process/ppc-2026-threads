#include "bortsova_a_integrals_rectangle/omp/include/ops_omp.hpp"

#include <cstdint>
#include <vector>

#include "bortsova_a_integrals_rectangle/common/include/common.hpp"

namespace bortsova_a_integrals_rectangle {

BortsovaAIntegralsRectangleOMP::BortsovaAIntegralsRectangleOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool BortsovaAIntegralsRectangleOMP::ValidationImpl() {
  const auto &input = GetInput();
  return input.func && !input.lower_bounds.empty() && input.lower_bounds.size() == input.upper_bounds.size() &&
         input.num_steps > 0;
}

bool BortsovaAIntegralsRectangleOMP::PreProcessingImpl() {
  const auto &input = GetInput();
  func_ = input.func;
  num_steps_ = input.num_steps;
  dims_ = static_cast<int>(input.lower_bounds.size());

  midpoints_.resize(dims_);
  volume_ = 1.0;
  total_points_ = 1;

  for (int di = 0; di < dims_; di++) {
    double step = (input.upper_bounds[di] - input.lower_bounds[di]) / static_cast<double>(num_steps_);
    volume_ *= step;
    total_points_ *= num_steps_;

    midpoints_[di].resize(num_steps_);
    for (int si = 0; si < num_steps_; si++) {
      midpoints_[di][si] = input.lower_bounds[di] + ((si + 0.5) * step);
    }
  }

  return true;
}

bool BortsovaAIntegralsRectangleOMP::RunImpl() {
  double sum = 0.0;
  const int dims = dims_;
  const int num_steps = num_steps_;
  const int64_t total_points = total_points_;
  const auto &midpoints = midpoints_;
  const auto &func = func_;

#pragma omp parallel default(none) shared(sum, dims, num_steps, total_points, midpoints, func)
  {
    std::vector<int> indices(dims, 0);
    std::vector<double> point(dims, 0.0);

#pragma omp for reduction(+ : sum)
    for (int64_t pt = 0; pt < total_points; pt++) {
      int64_t tmp = pt;
      for (int di = dims - 1; di >= 0; di--) {
        indices[di] = static_cast<int>(tmp % num_steps);
        tmp /= num_steps;
      }

      for (int di = 0; di < dims; di++) {
        point[di] = midpoints[di][indices[di]];
      }
      sum += func(point);
    }
  }

  GetOutput() = sum * volume_;
  return true;
}

bool BortsovaAIntegralsRectangleOMP::PostProcessingImpl() {
  return true;
}

}  // namespace bortsova_a_integrals_rectangle
