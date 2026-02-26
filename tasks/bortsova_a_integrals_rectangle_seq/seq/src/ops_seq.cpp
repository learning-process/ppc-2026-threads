#include "bortsova_a_integrals_rectangle_seq/seq/include/ops_seq.hpp"

#include <cstdint>
#include <vector>

#include "bortsova_a_integrals_rectangle_seq/common/include/common.hpp"

namespace bortsova_a_integrals_rectangle_seq {

BortsovaAIntegralsRectangleSEQ::BortsovaAIntegralsRectangleSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool BortsovaAIntegralsRectangleSEQ::ValidationImpl() {
  const auto &input = GetInput();
  return input.func && !input.lower_bounds.empty() && input.lower_bounds.size() == input.upper_bounds.size() &&
         input.num_steps > 0;
}

bool BortsovaAIntegralsRectangleSEQ::PreProcessingImpl() {
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

bool BortsovaAIntegralsRectangleSEQ::RunImpl() {
  double sum = 0.0;
  std::vector<int> indices(dims_, 0);
  std::vector<double> point(dims_);

  for (int64_t pt = 0; pt < total_points_; pt++) {
    for (int di = 0; di < dims_; di++) {
      point[di] = midpoints_[di][indices[di]];
    }
    sum += func_(point);

    for (int di = dims_ - 1; di >= 0; di--) {
      indices[di]++;
      if (indices[di] < num_steps_) {
        break;
      }
      indices[di] = 0;
    }
  }

  GetOutput() = sum * volume_;
  return true;
}

bool BortsovaAIntegralsRectangleSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace bortsova_a_integrals_rectangle_seq
