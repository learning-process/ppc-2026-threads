#include "bortsova_a_integrals_rectangle/tbb/include/ops_tbb.hpp"

#include <cstdint>
#include <functional>
#include <vector>

#include "bortsova_a_integrals_rectangle/common/include/common.hpp"
#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/parallel_reduce.h"

namespace bortsova_a_integrals_rectangle {

BortsovaAIntegralsRectangleTBB::BortsovaAIntegralsRectangleTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool BortsovaAIntegralsRectangleTBB::ValidationImpl() {
  const auto &input = GetInput();
  return input.func && !input.lower_bounds.empty() && input.lower_bounds.size() == input.upper_bounds.size() &&
         input.num_steps > 0;
}

bool BortsovaAIntegralsRectangleTBB::PreProcessingImpl() {
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

bool BortsovaAIntegralsRectangleTBB::RunImpl() {
  int dims = dims_;
  int num_steps = num_steps_;
  const auto &midpoints = midpoints_;
  const auto &func = func_;

  double sum = tbb::parallel_reduce(tbb::blocked_range<int64_t>(0, total_points_), 0.0,
                                    [&](const tbb::blocked_range<int64_t> &range, double local_sum) {
    std::vector<int> indices(dims, 0);
    std::vector<double> point(dims);

    int64_t temp = range.begin();
    for (int di = dims - 1; di >= 0; di--) {
      indices[di] = static_cast<int>(temp % num_steps);
      temp /= num_steps;
    }

    for (int64_t pt = range.begin(); pt < range.end(); pt++) {
      for (int di = 0; di < dims; di++) {
        point[di] = midpoints[di][indices[di]];
      }
      local_sum += func(point);

      for (int di = dims - 1; di >= 0; di--) {
        indices[di]++;
        if (indices[di] < num_steps) {
          break;
        }
        indices[di] = 0;
      }
    }
    return local_sum;
  }, std::plus<>());

  GetOutput() = sum * volume_;
  return true;
}

bool BortsovaAIntegralsRectangleTBB::PostProcessingImpl() {
  return true;
}

}  // namespace bortsova_a_integrals_rectangle
