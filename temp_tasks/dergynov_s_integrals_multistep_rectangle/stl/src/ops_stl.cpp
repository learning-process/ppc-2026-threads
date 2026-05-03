#include "dergynov_s_integrals_multistep_rectangle/stl/include/ops_stl.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <thread>
#include <utility>
#include <vector>

#include "dergynov_s_integrals_multistep_rectangle/common/include/common.hpp"

namespace dergynov_s_integrals_multistep_rectangle {
namespace {

bool ComputeSlice(const std::function<double(const std::vector<double> &)> &function,
                  const std::vector<std::pair<double, double>> &borders, const std::vector<double> &step_sizes,
                  size_t inner_dimensions, int64_t inner_cells, int steps_per_dimension, int outer_index,
                  std::vector<int> &indices, std::vector<double> &point, double &result) {
  indices[inner_dimensions] = outer_index;
  point[inner_dimensions] =
      borders[inner_dimensions].first + ((static_cast<double>(outer_index) + 0.5) * step_sizes[inner_dimensions]);

  for (int64_t cell = 0; cell < inner_cells; ++cell) {
    for (size_t dim_index = 0; dim_index < inner_dimensions; ++dim_index) {
      point[dim_index] =
          borders[dim_index].first + ((static_cast<double>(indices[dim_index]) + 0.5) * step_sizes[dim_index]);
    }

    double func_value = function(point);
    if (!std::isfinite(func_value)) {
      return false;
    }
    result += func_value;

    for (size_t position = 0; position < inner_dimensions; ++position) {
      if (++indices[position] < steps_per_dimension) {
        break;
      }
      indices[position] = 0;
    }
  }
  return true;
}

void ThreadWorker(const std::function<double(const std::vector<double> &)> &function,
                  const std::vector<std::pair<double, double>> &borders, const std::vector<double> &step_sizes,
                  size_t total_dimensions, size_t inner_dimensions, int64_t inner_cells, int steps_per_dimension,
                  int start_outer, int end_outer, double &partial_result, std::atomic<bool> &error_flag) {
  std::vector<int> indices(total_dimensions, 0);
  std::vector<double> point(total_dimensions);
  double local_sum = 0.0;

  for (int outer_idx = start_outer; outer_idx < end_outer && !error_flag.load(); ++outer_idx) {
    double slice_sum = 0.0;
    if (!ComputeSlice(function, borders, step_sizes, inner_dimensions, inner_cells, steps_per_dimension, outer_idx,
                      indices, point, slice_sum)) {
      error_flag.store(true);
      return;
    }
    local_sum += slice_sum;
  }

  partial_result = local_sum;
}

}  // namespace

DergynovSIntegralsMultistepRectangleSTL::DergynovSIntegralsMultistepRectangleSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0.0;
}

bool DergynovSIntegralsMultistepRectangleSTL::ValidationImpl() {
  const auto &[function, borders, steps_count] = GetInput();

  if (!function) {
    return false;
  }
  if (steps_count <= 0) {
    return false;
  }
  if (borders.empty()) {
    return false;
  }

  return std::ranges::all_of(borders, [](const auto &border) {
    return std::isfinite(border.first) && std::isfinite(border.second) && border.first < border.second;
  });
}

bool DergynovSIntegralsMultistepRectangleSTL::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

bool DergynovSIntegralsMultistepRectangleSTL::RunImpl() {
  const auto &[function, borders, steps_count] = GetInput();
  const size_t total_dimensions = borders.size();

  std::vector<double> step_sizes(total_dimensions);
  double cell_volume = 1.0;

  for (size_t dim_index = 0; dim_index < total_dimensions; ++dim_index) {
    step_sizes[dim_index] = (borders[dim_index].second - borders[dim_index].first) / static_cast<double>(steps_count);
    if (!(step_sizes[dim_index] > 0.0) || !std::isfinite(step_sizes[dim_index])) {
      return false;
    }
    cell_volume *= step_sizes[dim_index];
  }

  const size_t inner_dimensions = total_dimensions - 1;
  int64_t inner_cells = 1;
  for (size_t dim_index = 0; dim_index < inner_dimensions; ++dim_index) {
    inner_cells *= steps_count;
  }

  int hardware_threads = static_cast<int>(std::thread::hardware_concurrency());
  if (hardware_threads <= 0) {
    hardware_threads = 4;
  }
  const int actual_threads = std::min(hardware_threads, steps_count);

  std::vector<double> partial_results(actual_threads, 0.0);
  std::atomic<bool> error_flag{false};
  std::vector<std::thread> worker_threads;
  worker_threads.reserve(actual_threads);

  for (int thread_id = 0; thread_id < actual_threads; ++thread_id) {
    int chunk_size = steps_count / actual_threads;
    int remainder = steps_count % actual_threads;
    int start_index = (thread_id * chunk_size) + std::min(thread_id, remainder);
    int end_index = start_index + chunk_size + (thread_id < remainder ? 1 : 0);

    worker_threads.emplace_back([&, start_index, end_index, thread_id]() {
      ThreadWorker(function, borders, step_sizes, total_dimensions, inner_dimensions, inner_cells, steps_count,
                   start_index, end_index, partial_results[thread_id], error_flag);
    });
  }

  for (auto &worker : worker_threads) {
    worker.join();
  }

  if (error_flag.load()) {
    return false;
  }

  double total_sum = 0.0;
  for (double partial : partial_results) {
    total_sum += partial;
  }

  GetOutput() = total_sum * cell_volume;
  return std::isfinite(GetOutput());
}

bool DergynovSIntegralsMultistepRectangleSTL::PostProcessingImpl() {
  return std::isfinite(GetOutput());
}

}  // namespace dergynov_s_integrals_multistep_rectangle
