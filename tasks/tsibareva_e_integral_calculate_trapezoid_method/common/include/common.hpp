#pragma once

#include <cstdint>
#include <functional>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace tsibareva_e_integral_calculate_trapezoid_method {

enum class IntegralTestType : std::uint8_t {
  kSuccessSimple2D,
  kSuccessConstant2D,
  kSuccessSimple3D,
  kSuccessConstant3D,
  kInvalidLowerBoundEqual,
  kInvalidStepsNegative,
  kInvalidEmptyBounds,
};

struct Integral {
  std::vector<double> lo;
  std::vector<double> hi;
  std::vector<int> steps;
  std::function<double(const std::vector<double> &)> f;
  int dim{0};
};

using InType = Integral;
using OutType = double;
using TestType = std::tuple<IntegralTestType, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

inline Integral GenerateIntegralInput(IntegralTestType type) {
  Integral input;

  switch (type) {
    case IntegralTestType::kSuccessSimple2D: {
      input.dim = 2;
      input.lo = {0.0, 0.0};
      input.hi = {1.0, 1.0};
      input.steps = {100, 100};
      input.f = [](const std::vector<double> &x) { return (x[0] * x[0]) + (x[1] * x[1]); };  // x^2 + y^2
      break;
    }
    case IntegralTestType::kSuccessConstant2D: {
      input.dim = 2;
      input.lo = {0.0, 0.0};
      input.hi = {2.0, 3.0};
      input.steps = {50, 50};
      input.f = [](const std::vector<double> &) { return 5.0; };  // const
      break;
    }
    case IntegralTestType::kSuccessSimple3D: {
      input.dim = 3;
      input.lo = {0.0, 0.0, 0.0};
      input.hi = {1.0, 1.0, 1.0};
      input.steps = {50, 50, 50};
      input.f = [](const std::vector<double> &x) { return x[0] + x[1] + x[2]; };  // x + y + z
      break;
    }
    case IntegralTestType::kSuccessConstant3D: {
      input.dim = 3;
      input.lo = {0.0, 0.0, 0.0};
      input.hi = {2.0, 2.0, 2.0};
      input.steps = {40, 40, 40};
      input.f = [](const std::vector<double> &) { return 3.0; };
      break;
    }
    case IntegralTestType::kInvalidLowerBoundEqual: {
      input.dim = 2;
      input.lo = {1.0, 0.0};
      input.hi = {1.0, 1.0};
      input.steps = {10, 10};
      input.f = [](const std::vector<double> &x) { return x[0]; };
      break;
    }
    case IntegralTestType::kInvalidStepsNegative: {
      input.dim = 2;
      input.lo = {0.0, 0.0};
      input.hi = {1.0, 1.0};
      input.steps = {-5, 10};
      input.f = [](const std::vector<double> &x) { return x[0]; };
      break;
    }
    case IntegralTestType::kInvalidEmptyBounds: {
      input.dim = 0;
      input.lo = {};
      input.hi = {};
      input.steps = {};
      input.f = [](const std::vector<double> &) { return 0.0; };
      break;
    }
  }

  return input;
}

inline double GenerateExpectedOutput(IntegralTestType type) {
  switch (type) {
    case IntegralTestType::kSuccessSimple2D:
      return 2.0 / 3.0;
    case IntegralTestType::kSuccessConstant2D:
      return 30.0;
    case IntegralTestType::kSuccessSimple3D:
      return 1.5;
    case IntegralTestType::kSuccessConstant3D:
      return 24.0;
    case IntegralTestType::kInvalidLowerBoundEqual:
    case IntegralTestType::kInvalidStepsNegative:
    case IntegralTestType::kInvalidEmptyBounds:
      return 0.0;
  }
  return 0.0;
}

}  // namespace tsibareva_e_integral_calculate_trapezoid_method
