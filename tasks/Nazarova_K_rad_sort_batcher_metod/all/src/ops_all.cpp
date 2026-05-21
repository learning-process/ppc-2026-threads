#include "Nazarova_K_rad_sort_batcher_metod/all/include/ops_all.hpp"

#include <mpi.h>

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

NazarovaKCalcIntegRectanglesALL::NazarovaKCalcIntegRectanglesALL(const InType& in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool NazarovaKCalcIntegRectanglesALL::ValidationImpl() {
  return HasValidInput(GetInput());
}

bool NazarovaKCalcIntegRectanglesALL::PreProcessingImpl() {
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

bool NazarovaKCalcIntegRectanglesALL::RunImpl() {
  const auto& input = GetInput();

  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const std::size_t mpi_size = size > 0 ? static_cast<std::size_t>(size) : 1U;
  const std::size_t mpi_rank = rank >= 0 ? static_cast<std::size_t>(rank) : 0U;
  const std::size_t base_chunk = total_cells_ / mpi_size;
  const std::size_t remainder = total_cells_ % mpi_size;
  const std::size_t local_begin = (mpi_rank * base_chunk) + (mpi_rank < remainder ? mpi_rank : remainder);
  const std::size_t local_size = base_chunk + (mpi_rank < remainder ? 1U : 0U);
  const std::size_t local_end = local_begin + local_size;

  const double local_sum = tbb::parallel_reduce(
      tbb::blocked_range<std::size_t>(local_begin, local_end), 0.0,
      [&](const tbb::blocked_range<std::size_t>& range, double partial_sum) {
        std::vector<double> point(dimension_, 0.0);
        for (std::size_t linear_index = range.begin(); linear_index < range.end(); ++linear_index) {
          FillPointFromLinearIndex(input, step_sizes_, linear_index, point);
          partial_sum += input.function(point);
        }
        return partial_sum;
      },
      std::plus<>());

  double global_sum = 0.0;
  MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  result_ = global_sum * cell_volume_;
  GetOutput() = result_;
  return true;
}

bool NazarovaKCalcIntegRectanglesALL::PostProcessingImpl() {
  return true;
}

}  // namespace nazarova_k_calc_integ_rectangles
