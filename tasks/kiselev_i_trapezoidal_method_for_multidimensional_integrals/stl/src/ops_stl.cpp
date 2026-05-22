#include "kiselev_i_trapezoidal_method_for_multidimensional_integrals/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <future>
#include <vector>

#include "kiselev_i_trapezoidal_method_for_multidimensional_integrals/common/include/common.hpp"

namespace kiselev_i_trapezoidal_method_for_multidimensional_integrals {

namespace {

double ComputeChunk(const InType &input_data, const std::vector<int> &steps, int start_index, int end_index) {
  const double hx =
      static_cast<double>(input_data.right_bounds[0] - input_data.left_bounds[0]) / static_cast<double>(steps[0]);

  const double hy =
      static_cast<double>(input_data.right_bounds[1] - input_data.left_bounds[1]) / static_cast<double>(steps[1]);

  double local_result = 0.0;

  for (int x_index = start_index; x_index < end_index; x_index++) {
    const double x = input_data.left_bounds[0] + (static_cast<double>(x_index) * hx);

    const double weight_x = (x_index == 0 || x_index == steps[0]) ? 0.5 : 1.0;

    for (int y_index = 0; y_index <= steps[1]; y_index++) {
      const double y = input_data.left_bounds[1] + (static_cast<double>(y_index) * hy);

      const double weight_y = (y_index == 0 || y_index == steps[1]) ? 0.5 : 1.0;

      double value = 0.0;
      value = KiselevITestTaskSTL::FunctionTypeChoose(input_data.type_function, x, y);

      local_result += weight_x * weight_y * value;
    }
  }

  return local_result;
}

}  // namespace

KiselevITestTaskSTL::KiselevITestTaskSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());

  GetInput() = in;

  GetOutput() = 0.0;
}

bool KiselevITestTaskSTL::ValidationImpl() {
  return true;
}

bool KiselevITestTaskSTL::PreProcessingImpl() {
  GetOutput() = 0.0;

  return true;
}

double KiselevITestTaskSTL::FunctionTypeChoose(int type_x, double x, double y) {
  switch (type_x) {
    case 0:
      return (x * x) + (y * y);

    case 1:
      return std::sin(x) * std::cos(y);

    case 2:
      return std::sin(x) + std::cos(y);

    case 3:
      return std::exp(x + y);

    default:
      return x + y;
  }
}

double KiselevITestTaskSTL::ComputeIntegral(const std::vector<int> &steps) {
  const auto &input_data = GetInput();

  const int total_iterations = steps[0] + 1;

  int num_threads = 4;

  num_threads = std::min(num_threads, total_iterations);

  const int chunk_size = total_iterations / num_threads;

  const int remainder = total_iterations % num_threads;

  std::vector<std::future<double>> futures;

  int start_index = 0;

  for (int thread_index = 0; thread_index < num_threads; thread_index++) {
    int end_index = start_index + chunk_size;

    if (thread_index < remainder) {
      end_index++;
    }

    futures.emplace_back(
        std::async(std::launch::async, ComputeChunk, std::cref(input_data), std::cref(steps), start_index, end_index));

    start_index = end_index;
  }

  const double hx =
      static_cast<double>(input_data.right_bounds[0] - input_data.left_bounds[0]) / static_cast<double>(steps[0]);

  const double hy =
      static_cast<double>(input_data.right_bounds[1] - input_data.left_bounds[1]) / static_cast<double>(steps[1]);

  double result = 0.0;

  for (auto &future_result : futures) {
    result += future_result.get();
  }

  return result * hx * hy;
}

bool KiselevITestTaskSTL::RunImpl() {
  const auto &input_data = GetInput();

  if (input_data.left_bounds.size() != 2 || input_data.right_bounds.size() != 2 || input_data.step_n_size.size() != 2) {
    GetOutput() = 0.0;

    return true;
  }

  std::vector<int> steps = input_data.step_n_size;

  for (const auto &step_value : steps) {
    if (step_value <= 0) {
      GetOutput() = 0.0;

      return true;
    }
  }

  const double epsilon = input_data.epsilon;

  if (epsilon <= 0.0) {
    GetOutput() = ComputeIntegral(steps);

    return true;
  }

  double previous_result = ComputeIntegral(steps);

  double current_result = previous_result;

  const int max_iterations = 1;

  for (int iteration_index = 0; iteration_index < max_iterations; iteration_index++) {
    for (auto &step_value : steps) {
      step_value *= 2;
    }

    current_result = ComputeIntegral(steps);

    if (std::abs(current_result - previous_result) < epsilon) {
      break;
    }

    previous_result = current_result;
  }

  GetOutput() = current_result;

  return true;
}

bool KiselevITestTaskSTL::PostProcessingImpl() {
  return true;
}

}  // namespace kiselev_i_trapezoidal_method_for_multidimensional_integrals
