#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace shilin_n_monte_carlo_integration {

enum class FuncType : std::uint8_t {
  kConstant = 0,
  kLinear = 1,
  kProduct = 2,
  kSumSquares = 3,
  kSinProduct = 4,
};

class IntegrandFunction {
 public:
  static double Evaluate(FuncType func_type, const std::vector<double> &point) {
    switch (func_type) {
      case FuncType::kConstant:
        return 1.0;
      case FuncType::kLinear: {
        double sum = 0.0;
        for (double x : point) {
          sum += x;
        }
        return sum;
      }
      case FuncType::kProduct: {
        double prod = 1.0;
        for (double x : point) {
          prod *= x;
        }
        return prod;
      }
      case FuncType::kSumSquares: {
        double sum = 0.0;
        for (double x : point) {
          sum += x * x;
        }
        return sum;
      }
      case FuncType::kSinProduct: {
        double prod = 1.0;
        for (double x : point) {
          prod *= std::sin(x);
        }
        return prod;
      }
      default:
        return 0.0;
    }
  }

  static double AnalyticalIntegral(FuncType func_type, const std::vector<double> &lower,
                                   const std::vector<double> &upper) {
    size_t dim = lower.size();
    switch (func_type) {
      case FuncType::kConstant:
        return ComputeVolume(lower, upper);
      case FuncType::kLinear:
        return ComputeLinearIntegral(lower, upper, dim);
      case FuncType::kProduct:
        return ComputeProductIntegral(lower, upper, dim);
      case FuncType::kSumSquares:
        return ComputeSumSquaresIntegral(lower, upper, dim);
      case FuncType::kSinProduct:
        return ComputeSinProductIntegral(lower, upper, dim);
      default:
        return 0.0;
    }
  }

 private:
  static double ComputeVolume(const std::vector<double> &lower, const std::vector<double> &upper) {
    double vol = 1.0;
    for (size_t i = 0; i < lower.size(); ++i) {
      vol *= (upper[i] - lower[i]);
    }
    return vol;
  }

  static double VolumeExcludingDim(const std::vector<double> &lower, const std::vector<double> &upper, size_t exclude) {
    double product = 1.0;
    for (size_t j = 0; j < lower.size(); ++j) {
      if (j != exclude) {
        product *= (upper[j] - lower[j]);
      }
    }
    return product;
  }

  static double ComputeLinearIntegral(const std::vector<double> &lower, const std::vector<double> &upper, size_t dim) {
    double result = 0.0;
    for (size_t i = 0; i < dim; ++i) {
      double term = (upper[i] * upper[i] - lower[i] * lower[i]) / 2.0;
      result += term * VolumeExcludingDim(lower, upper, i);
    }
    return result;
  }

  static double ComputeProductIntegral(const std::vector<double> &lower, const std::vector<double> &upper, size_t dim) {
    double prod = 1.0;
    for (size_t i = 0; i < dim; ++i) {
      prod *= (upper[i] * upper[i] - lower[i] * lower[i]) / 2.0;
    }
    return prod;
  }

  static double ComputeSumSquaresIntegral(const std::vector<double> &lower, const std::vector<double> &upper,
                                          size_t dim) {
    double result = 0.0;
    for (size_t i = 0; i < dim; ++i) {
      double term = (upper[i] * upper[i] * upper[i] - lower[i] * lower[i] * lower[i]) / 3.0;
      result += term * VolumeExcludingDim(lower, upper, i);
    }
    return result;
  }

  static double ComputeSinProductIntegral(const std::vector<double> &lower, const std::vector<double> &upper,
                                          size_t dim) {
    double prod = 1.0;
    for (size_t i = 0; i < dim; ++i) {
      prod *= (-std::cos(upper[i]) + std::cos(lower[i]));
    }
    return prod;
  }
};

using InType = std::tuple<std::vector<double>, std::vector<double>, int, FuncType>;
using OutType = double;
using TestType = std::tuple<InType, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace shilin_n_monte_carlo_integration
