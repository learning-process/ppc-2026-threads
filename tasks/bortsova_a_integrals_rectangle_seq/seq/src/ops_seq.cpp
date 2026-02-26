#include "bortsova_a_integrals_rectangle_seq/seq/include/ops_seq.hpp"

#include <cstdint>
#include <vector>

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
  lower_bounds_ = input.lower_bounds;
  upper_bounds_ = input.upper_bounds;
  num_steps_ = input.num_steps;
  return true;
}

bool BortsovaAIntegralsRectangleSEQ::RunImpl() {
  int dims = static_cast<int>(lower_bounds_.size());

  std::vector<double> step_sizes(dims);
  for (int i = 0; i < dims; i++) {
    step_sizes[i] = (upper_bounds_[i] - lower_bounds_[i]) / static_cast<double>(num_steps_);
  }

  int64_t total_points = 1;
  for (int i = 0; i < dims; i++) {
    total_points *= num_steps_;
  }

  double sum = 0.0;
  std::vector<int> indices(dims, 0);
  std::vector<double> point(dims);

  for (int64_t p = 0; p < total_points; p++) {
    for (int d = 0; d < dims; d++) {
      point[d] = lower_bounds_[d] + (indices[d] + 0.5) * step_sizes[d];
    }
    sum += func_(point);

    for (int d = dims - 1; d >= 0; d--) {
      indices[d]++;
      if (indices[d] < num_steps_) {
        break;
      }
      indices[d] = 0;
    }
  }

  double volume = 1.0;
  for (int i = 0; i < dims; i++) {
    volume *= step_sizes[i];
  }

  GetOutput() = sum * volume;
  return true;
}

bool BortsovaAIntegralsRectangleSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace bortsova_a_integrals_rectangle_seq
