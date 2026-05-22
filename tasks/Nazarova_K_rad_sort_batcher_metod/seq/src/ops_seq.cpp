#include "Nazarova_K_rad_sort_batcher_metod/seq/include/ops_seq.hpp"

#include <cmath>
#include <cstddef>
#include <vector>

#include "Nazarova_K_rad_sort_batcher_metod/common/include/common.hpp"

namespace nazarova_k_calc_integ_rectangles {

NazarovaKCalcIntegRectanglesSEQ::NazarovaKCalcIntegRectanglesSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool NazarovaKCalcIntegRectanglesSEQ::HasValidInput() {
  const auto &input = GetInput();
  const std::size_t dimension = input.lower_bounds.size();
  if (!input.function || dimension == 0U || input.upper_bounds.size() != dimension || input.steps.size() != dimension) {
    return false;
  }

  for (std::size_t i = 0; i < dimension; ++i) {
    if (input.steps[i] == 0U || !std::isfinite(input.lower_bounds[i]) || !std::isfinite(input.upper_bounds[i]) ||
        input.lower_bounds[i] > input.upper_bounds[i]) {
      return false;
    }
  }

  return true;
}

bool NazarovaKCalcIntegRectanglesSEQ::ValidationImpl() {
  return HasValidInput();
}

bool NazarovaKCalcIntegRectanglesSEQ::PreProcessingImpl() {
  const auto &input = GetInput();
  dimension_ = input.lower_bounds.size();
  step_sizes_.assign(dimension_, 0.0);
  cell_volume_ = 1.0;
  result_ = 0.0;

  for (std::size_t i = 0; i < dimension_; ++i) {
    step_sizes_[i] = (input.upper_bounds[i] - input.lower_bounds[i]) / static_cast<double>(input.steps[i]);
    cell_volume_ *= step_sizes_[i];
  }

  GetOutput() = 0.0;
  return true;
}

bool NazarovaKCalcIntegRectanglesSEQ::RunImpl() {
  const auto &input = GetInput();
  std::vector<double> point(dimension_, 0.0);
  std::vector<std::size_t> indices(dimension_, 0U);
  double sum = 0.0;

  while (true) {
    for (std::size_t i = 0; i < dimension_; ++i) {
      point[i] = input.lower_bounds[i] + ((static_cast<double>(indices[i]) + 0.5) * step_sizes_[i]);
    }
    sum += input.function(point);

    std::size_t axis = 0U;
    while (axis < dimension_) {
      ++indices[axis];
      if (indices[axis] < input.steps[axis]) {
        break;
      }
      indices[axis] = 0U;
      ++axis;
    }
    if (axis == dimension_) {
      break;
    }
  }

  result_ = sum * cell_volume_;
  GetOutput() = result_;
  return true;
}

bool NazarovaKCalcIntegRectanglesSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace nazarova_k_calc_integ_rectangles
