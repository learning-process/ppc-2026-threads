#include "Nazarova_K_rad_sort_batcher_metod/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <future>
#include <limits>
#include <thread>
#include <vector>

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

NazarovaKCalcIntegRectanglesSTL::NazarovaKCalcIntegRectanglesSTL(const InType& in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool NazarovaKCalcIntegRectanglesSTL::ValidationImpl() {
  return HasValidInput(GetInput());
}

bool NazarovaKCalcIntegRectanglesSTL::PreProcessingImpl() {
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

bool NazarovaKCalcIntegRectanglesSTL::RunImpl() {
  const auto& input = GetInput();
  std::size_t thread_count = std::thread::hardware_concurrency();
  if (thread_count == 0U) {
    thread_count = 1U;
  }
  thread_count = std::min(thread_count, total_cells_ == 0U ? 1U : total_cells_);

  const std::size_t base_chunk = total_cells_ / thread_count;
  const std::size_t remainder = total_cells_ % thread_count;

  std::vector<std::future<double>> futures;
  futures.reserve(thread_count);

  std::size_t begin = 0U;
  for (std::size_t thread_id = 0U; thread_id < thread_count; ++thread_id) {
    const std::size_t extra = thread_id < remainder ? 1U : 0U;
    const std::size_t end = begin + base_chunk + extra;
    futures.emplace_back(std::async(std::launch::async, [&, begin, end] {
      std::vector<double> point(dimension_, 0.0);
      double local_sum = 0.0;
      for (std::size_t linear_index = begin; linear_index < end; ++linear_index) {
        FillPointFromLinearIndex(input, step_sizes_, linear_index, point);
        local_sum += input.function(point);
      }
      return local_sum;
    }));
    begin = end;
  }

  double sum = 0.0;
  for (auto& future : futures) {
    sum += future.get();
  }

  result_ = sum * cell_volume_;
  GetOutput() = result_;
  return true;
}

bool NazarovaKCalcIntegRectanglesSTL::PostProcessingImpl() {
  return true;
}

}  // namespace nazarova_k_calc_integ_rectangles
