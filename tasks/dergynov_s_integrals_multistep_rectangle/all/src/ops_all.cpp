#include "dergynov_s_integrals_multistep_rectangle/all/include/ops_all.hpp"

#include <omp.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <functional>
#include <thread>
#include <utility>
#include <vector>

#include "dergynov_s_integrals_multistep_rectangle/common/include/common.hpp"

namespace dergynov_s_integrals_multistep_rectangle {
namespace {

bool ValidateBorders(const std::vector<std::pair<double, double>> &borders) {
  return std::ranges::all_of(borders, [](const auto &border) {
    return std::isfinite(border.first) && std::isfinite(border.second) && border.first < border.second;
  });
}

void ComputePoint(size_t linear_idx, int steps_count, int dimensions,
                  const std::vector<std::pair<double, double>> &limits, const std::vector<double> &step_sizes,
                  std::vector<double> &coordinates) {
  size_t temp = linear_idx;
  for (int dimension = dimensions - 1; dimension >= 0; --dimension) {
    int idx_val = static_cast<int>(temp % static_cast<size_t>(steps_count));
    temp /= static_cast<size_t>(steps_count);
    coordinates[dimension] = limits[dimension].first + ((static_cast<double>(idx_val) + 0.5) * step_sizes[dimension]);
  }
}

double ComputeSEQ(const std::function<double(const std::vector<double> &)> &function,
                  const std::vector<std::pair<double, double>> &limits, int steps_count, int dimensions,
                  const std::vector<double> &step_sizes) {
  size_t total_points = 1;
  for (int dim_idx = 0; dim_idx < dimensions; ++dim_idx) {
    total_points *= static_cast<size_t>(steps_count);
  }

  double total_sum = 0.0;
  for (size_t linear_index = 0; linear_index < total_points; ++linear_index) {
    std::vector<double> point(dimensions);
    ComputePoint(linear_index, steps_count, dimensions, limits, step_sizes, point);
    total_sum += function(point);
  }
  return total_sum;
}

double ComputeOMP(const std::function<double(const std::vector<double> &)> &function,
                  const std::vector<std::pair<double, double>> &limits, int steps_count, int dimensions,
                  const std::vector<double> &step_sizes) {
  size_t total_points = 1;
  for (int dim_idx = 0; dim_idx < dimensions; ++dim_idx) {
    total_points *= static_cast<size_t>(steps_count);
  }

  std::vector<double> thread_sums(omp_get_max_threads(), 0.0);
  int error_flag_value = 0;

#pragma omp parallel default(none) \
    shared(function, limits, step_sizes, steps_count, dimensions, total_points, thread_sums, error_flag_value)
  {
    int thread_id = omp_get_thread_num();
    double local_sum = 0.0;

#pragma omp for schedule(static)
    for (size_t linear_index = 0; linear_index < total_points; ++linear_index) {
      if (error_flag_value != 0) {
        continue;
      }

      std::vector<double> point(dimensions);
      ComputePoint(linear_index, steps_count, dimensions, limits, step_sizes, point);

      double func_value = function(point);
      if (!std::isfinite(func_value)) {
#pragma omp atomic write
        error_flag_value = 1;
        continue;
      }
      local_sum += func_value;
    }

    thread_sums[thread_id] = local_sum;
  }

  if (error_flag_value != 0) {
    return 0.0;
  }

  double total_sum = 0.0;
  for (double sum_value : thread_sums) {
    total_sum += sum_value;
  }
  return total_sum;
}

double ComputeTBB(const std::function<double(const std::vector<double> &)> &function,
                  const std::vector<std::pair<double, double>> &limits, int steps_count, int dimensions,
                  const std::vector<double> &step_sizes) {
  size_t total_points = 1;
  for (int dim_idx = 0; dim_idx < dimensions; ++dim_idx) {
    total_points *= static_cast<size_t>(steps_count);
  }

  double total_sum = tbb::parallel_reduce(tbb::blocked_range<size_t>(0, total_points), 0.0,
                                          [&](const tbb::blocked_range<size_t> &range, double local_sum) {
    for (size_t linear_index = range.begin(); linear_index != range.end(); ++linear_index) {
      std::vector<double> point(dimensions);
      ComputePoint(linear_index, steps_count, dimensions, limits, step_sizes, point);
      local_sum += function(point);
    }
    return local_sum;
  }, [](double first, double second) { return first + second; });

  return total_sum;
}

void ProcessRangeSTL(size_t start, size_t end, const std::function<double(const std::vector<double> &)> &function,
                     const std::vector<std::pair<double, double>> &limits, int steps_count, int dimensions,
                     const std::vector<double> &step_sizes, std::atomic<bool> &error_occurred, double &result) {
  double local_sum = 0.0;
  for (size_t linear_index = start; linear_index < end && !error_occurred.load(); ++linear_index) {
    std::vector<double> point(dimensions);
    ComputePoint(linear_index, steps_count, dimensions, limits, step_sizes, point);
    double func_value = function(point);
    if (!std::isfinite(func_value)) {
      error_occurred.store(true);
      return;
    }
    local_sum += func_value;
  }
  result = local_sum;
}

double ComputeSTL(const std::function<double(const std::vector<double> &)> &function,
                  const std::vector<std::pair<double, double>> &limits, int steps_count, int dimensions,
                  const std::vector<double> &step_sizes) {
  size_t total_points = 1;
  for (int dim_idx = 0; dim_idx < dimensions; ++dim_idx) {
    total_points *= static_cast<size_t>(steps_count);
  }

  unsigned int thread_count = std::thread::hardware_concurrency();
  if (thread_count == 0) {
    thread_count = 1;
  }
  thread_count = std::min(thread_count, 2U);
  thread_count = std::min(thread_count, static_cast<unsigned int>(total_points));
  if (thread_count == 0) {
    thread_count = 1;
  }

  std::vector<double> partial_results(thread_count, 0.0);
  std::atomic<bool> error_occurred{false};
  std::vector<std::thread> worker_threads;
  worker_threads.reserve(thread_count);

  size_t chunk_size = total_points / thread_count;
  size_t remainder = total_points % thread_count;
  size_t current_start = 0;

  for (unsigned int thread_index = 0; thread_index < thread_count; ++thread_index) {
    size_t current_end = current_start + chunk_size + (thread_index < remainder ? 1 : 0);
    worker_threads.emplace_back(ProcessRangeSTL, current_start, current_end, std::cref(function), std::cref(limits),
                                steps_count, dimensions, std::cref(step_sizes), std::ref(error_occurred),
                                std::ref(partial_results[thread_index]));
    current_start = current_end;
  }

  for (auto &worker : worker_threads) {
    worker.join();
  }

  if (error_occurred.load()) {
    return 0.0;
  }

  double total_sum = 0.0;
  for (double partial_value : partial_results) {
    total_sum += partial_value;
  }
  return total_sum;
}

double ComputeALL(const std::function<double(const std::vector<double> &)> &function,
                  const std::vector<std::pair<double, double>> &limits, int steps_count, int dimensions) {
  std::vector<double> step_sizes(dimensions);
  double cell_volume = 1.0;

  for (int dim_idx = 0; dim_idx < dimensions; ++dim_idx) {
    step_sizes[dim_idx] = (limits[dim_idx].second - limits[dim_idx].first) / static_cast<double>(steps_count);
    cell_volume *= step_sizes[dim_idx];
  }

  double seq_result = ComputeSEQ(function, limits, steps_count, dimensions, step_sizes);
  double omp_result = ComputeOMP(function, limits, steps_count, dimensions, step_sizes);
  double stl_result = ComputeSTL(function, limits, steps_count, dimensions, step_sizes);
  double tbb_result = ComputeTBB(function, limits, steps_count, dimensions, step_sizes);

  const double epsilon = 1e-6;
  bool results_match = true;

  if (std::abs(seq_result - omp_result) > epsilon) {
    results_match = false;
  }
  if (std::abs(seq_result - stl_result) > epsilon) {
    results_match = false;
  }
  if (std::abs(seq_result - tbb_result) > epsilon) {
    results_match = false;
  }

  if (!results_match) {
    return 0.0;
  }

  return tbb_result * cell_volume;
}

}  // namespace

DergynovSIntegralsMultistepRectangleALL::DergynovSIntegralsMultistepRectangleALL(const InType &input_params) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = input_params;
  GetOutput() = 0.0;
}

bool DergynovSIntegralsMultistepRectangleALL::ValidationImpl() {
  const auto &[target_function, integration_limits, step_count] = GetInput();

  if (!target_function) {
    return false;
  }
  if (step_count <= 0) {
    return false;
  }
  if (integration_limits.empty()) {
    return false;
  }

  return ValidateBorders(integration_limits);
}

bool DergynovSIntegralsMultistepRectangleALL::PreProcessingImpl() {
  GetOutput() = 0.0;
  return true;
}

bool DergynovSIntegralsMultistepRectangleALL::RunImpl() {
  const auto &input_data = GetInput();
  const auto &target_function = std::get<0>(input_data);
  const auto &integration_limits = std::get<1>(input_data);
  const int step_count = std::get<2>(input_data);
  const int total_dimensions = static_cast<int>(integration_limits.size());

  double final_result = ComputeALL(target_function, integration_limits, step_count, total_dimensions);

  if (final_result == 0.0) {
    return false;
  }

  GetOutput() = final_result;
  return std::isfinite(GetOutput());
}

bool DergynovSIntegralsMultistepRectangleALL::PostProcessingImpl() {
  return std::isfinite(GetOutput());
}

}  // namespace dergynov_s_integrals_multistep_rectangle
