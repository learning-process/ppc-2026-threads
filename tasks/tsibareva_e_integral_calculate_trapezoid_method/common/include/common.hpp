#pragma once

#include <cmath>
#include <functional>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace tsibareva_e_integral_calculate_trapezoid_method {

enum class IntegralTestType {
  SUCCESS_SIMPLE_2D,
  SUCCESS_CONSTANT_2D,
  SUCCESS_SIMPLE_3D,
  SUCCESS_CONSTANT_3D,
  INVALID_LOWER_BOUND_EQUAL,
  INVALID_STEPS_NEGATIVE,
  INVALID_EMPTY_BOUNDS
};

struct IntegralInput {
  std::vector<double> lower_bounds;
  std::vector<double> upper_bounds;
  std::vector<int> num_steps;
  std::function<double(const std::vector<double>&)> function;
  int dimension;
};

using InType = IntegralInput;
using OutType = double;
using TestType = std::tuple<IntegralTestType, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

inline IntegralInput GenerateIntegralInput(IntegralTestType type) {
  IntegralInput input;
  
  switch (type) {
    case IntegralTestType::SUCCESS_SIMPLE_2D: {
      input.dimension = 2;
      input.lower_bounds = {0.0, 0.0};
      input.upper_bounds = {1.0, 1.0};
      input.num_steps = {100, 100};
      input.function = [](const std::vector<double>& x) { return x[0] * x[0] + x[1] * x[1]; };
      break;
    }
    case IntegralTestType::SUCCESS_CONSTANT_2D: {
      input.dimension = 2;
      input.lower_bounds = {0.0, 0.0};
      input.upper_bounds = {2.0, 3.0};
      input.num_steps = {50, 50};
      input.function = [](const std::vector<double>&) { return 5.0; };
      break;
    }
    case IntegralTestType::SUCCESS_SIMPLE_3D: {
      input.dimension = 3;
      input.lower_bounds = {0.0, 0.0, 0.0};
      input.upper_bounds = {1.0, 1.0, 1.0};
      input.num_steps = {50, 50, 50};
      input.function = [](const std::vector<double>& x) { return x[0] + x[1] + x[2]; };
      break;
    }
    case IntegralTestType::SUCCESS_CONSTANT_3D: {
      input.dimension = 3;
      input.lower_bounds = {0.0, 0.0, 0.0};
      input.upper_bounds = {2.0, 2.0, 2.0};
      input.num_steps = {40, 40, 40};
      input.function = [](const std::vector<double>&) { return 3.0; };
      break;
    }
    case IntegralTestType::INVALID_LOWER_BOUND_EQUAL: {
      input.dimension = 2;
      input.lower_bounds = {1.0, 0.0};
      input.upper_bounds = {1.0, 1.0};
      input.num_steps = {10, 10};
      input.function = [](const std::vector<double>& x) { return x[0]; };
      break;
    }
    case IntegralTestType::INVALID_STEPS_NEGATIVE: {
      input.dimension = 2;
      input.lower_bounds = {0.0, 0.0};
      input.upper_bounds = {1.0, 1.0};
      input.num_steps = {-5, 10};
      input.function = [](const std::vector<double>& x) { return x[0]; };
      break;
    }
    case IntegralTestType::INVALID_EMPTY_BOUNDS: {
      input.dimension = 0;
      input.lower_bounds = {};
      input.upper_bounds = {};
      input.num_steps = {};
      input.function = [](const std::vector<double>&) { return 0.0; };
      break;
    }
  }
  
  return input;
}

inline double GenerateExpectedOutput(IntegralTestType type) {
  switch (type) {
    case IntegralTestType::SUCCESS_SIMPLE_2D:
      return 2.0 / 3.0;
    case IntegralTestType::SUCCESS_CONSTANT_2D:
      return 30.0;
    case IntegralTestType::SUCCESS_SIMPLE_3D:
      return 1.5;
    case IntegralTestType::SUCCESS_CONSTANT_3D:
      return 24.0;
    case IntegralTestType::INVALID_LOWER_BOUND_EQUAL:
    case IntegralTestType::INVALID_STEPS_NEGATIVE:
    case IntegralTestType::INVALID_EMPTY_BOUNDS:
      return 0.0;
  }
  return 0.0;
}

}  // namespace tsibareva_e_integral_calculate_trapezoid_method