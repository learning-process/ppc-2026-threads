#include "Nazarova_K_rad_sort_batcher_metod/tbb/include/ops_tbb.hpp"

#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <vector>

#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>

#include "Nazarova_K_rad_sort_batcher_metod/common/include/common.hpp"

namespace nazarova_k_calc_integ_rectangles {

namespace {

bool HasValidInput(const InType& input) {
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

void FillPointFromLinearIndex(const InType& input, const std::vector<double>& step_sizes, std::size_t linear_index,
                              std::vector<double>& point) {
  std::size_t current_index = linear_index;
  for (std::size_t axis = 0U; axis < point.size(); ++axis) {
    const std::size_t coordinate_index = current_index % input.steps[axis];
    current_index /= input.steps[axis];
    point[axis] = input.lower_bounds[axis] + ((static_cast<double>(coordinate_index) + 0.5) * step_sizes[axis]);
  }
}

}  // namespace

NazarovaKCalcIntegRectanglesTBB::NazarovaKCalcIntegRectanglesTBB(const InType& in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool NazarovaKCalcIntegRectanglesTBB::ValidationImpl() {
  return HasValidInput(GetInput());
}

bool NazarovaKCalcIntegRectanglesTBB::PreProcessingImpl() {
  const auto& input = GetInput();
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

bool NazarovaKCalcIntegRectanglesTBB::RunImpl() {
  const auto& input = GetInput();
  const double sum = tbb::parallel_reduce(
      tbb::blocked_range<std::size_t>(0U, total_cells_), 0.0,
      [&](const tbb::blocked_range<std::size_t>& range, double local_sum) {
        std::vector<double> point(dimension_, 0.0);
        for (std::size_t linear_index = range.begin(); linear_index < range.end(); ++linear_index) {
          FillPointFromLinearIndex(input, step_sizes_, linear_index, point);
          local_sum += input.function(point);
        }
        return local_sum;
      },
      std::plus<>());

  result_ = sum * cell_volume_;
  GetOutput() = result_;
  return true;
}

bool NazarovaKCalcIntegRectanglesTBB::PostProcessingImpl() {
  return true;
}

}  // namespace nazarova_k_calc_integ_rectangles
