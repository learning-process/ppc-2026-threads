#include "Nazarova_K_rad_sort_batcher_metod/omp/include/ops_omp.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "Nazarova_K_rad_sort_batcher_metod/common/include/common.hpp"

namespace nazarova_k_calc_integ_rectangles {

NazarovaKCalcIntegRectanglesOMP::NazarovaKCalcIntegRectanglesOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool NazarovaKCalcIntegRectanglesOMP::HasValidInput() {
  const auto &input = GetInput();
  const std::size_t dimension = input.lower_bounds.size();
  if (!input.function || dimension == 0U || input.upper_bounds.size() != dimension || input.steps.size() != dimension) {
    return false;
  }

  std::size_t total_cells = 1U;
  for (std::size_t i = 0; i < dimension; ++i) {
    if (input.steps[i] == 0U || !std::isfinite(input.lower_bounds[i]) || !std::isfinite(input.upper_bounds[i]) ||
        input.lower_bounds[i] > input.upper_bounds[i]) {
      return false;
    }
    if (total_cells > (std::numeric_limits<std::size_t>::max() / input.steps[i])) {
      return false;
    }
    total_cells *= input.steps[i];
  }

  return true;
}

bool NazarovaKCalcIntegRectanglesOMP::ValidationImpl() {
  return HasValidInput();
}

bool NazarovaKCalcIntegRectanglesOMP::PreProcessingImpl() {
  const auto &input = GetInput();
  dimension_ = input.lower_bounds.size();
  step_sizes_.assign(dimension_, 0.0);
  cell_volume_ = 1.0;
  total_cells_ = 1U;
  result_ = 0.0;

  for (std::size_t i = 0; i < dimension_; ++i) {
    step_sizes_[i] = (input.upper_bounds[i] - input.lower_bounds[i]) / static_cast<double>(input.steps[i]);
    cell_volume_ *= step_sizes_[i];
    total_cells_ *= input.steps[i];
  }

  GetOutput() = 0.0;
  return true;
}

bool NazarovaKCalcIntegRectanglesOMP::RunImpl() {
  const auto &input = GetInput();
  const std::size_t dimension = dimension_;
  const std::size_t total_cells = total_cells_;
  const auto &step_sizes = step_sizes_;
  double sum = 0.0;

#pragma omp parallel default(none) shared(input, step_sizes) firstprivate(dimension, total_cells) reduction(+ : sum)
  {
    std::vector<double> point(dimension, 0.0);

#pragma omp for schedule(static)
    for (std::size_t linear_index = 0U; linear_index < total_cells; ++linear_index) {
      std::size_t current_index = linear_index;
      for (std::size_t axis = 0U; axis < dimension; ++axis) {
        const std::size_t coordinate_index = current_index % input.steps[axis];
        current_index /= input.steps[axis];
        point[axis] = input.lower_bounds[axis] + ((static_cast<double>(coordinate_index) + 0.5) * step_sizes[axis]);
      }
      sum += input.function(point);
    }
  }

  result_ = sum * cell_volume_;
  GetOutput() = result_;
  return true;
}

bool NazarovaKCalcIntegRectanglesOMP::PostProcessingImpl() {
  return true;
}

}  // namespace nazarova_k_calc_integ_rectangles
