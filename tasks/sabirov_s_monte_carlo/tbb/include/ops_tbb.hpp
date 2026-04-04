#pragma once

#include <cmath>
#include <vector>

#include "sabirov_s_monte_carlo/common/include/common.hpp"
#include "task/include/task.hpp"

namespace sabirov_s_monte_carlo {

namespace detail {

inline double EvalLinear(const std::vector<double> &point) {
  double s = 0.0;
  for (double x : point) {
    s += x;
  }
  return s;
}

inline double EvalSumCubes(const std::vector<double> &point) {
  double s = 0.0;
  for (double x : point) {
    s += x * x * x;
  }
  return s;
}

inline double EvalCosProduct(const std::vector<double> &point) {
  double p = 1.0;
  for (double x : point) {
    p *= std::cos(x);
  }
  return p;
}

inline double EvalExpNeg(const std::vector<double> &point) {
  double s = 0.0;
  for (double x : point) {
    s += x;
  }
  return std::exp(-s);
}

inline double EvalMixedPoly(const std::vector<double> &point) {
  double s = 0.0;
  for (double x : point) {
    s += (x * x) + x;
  }
  return s;
}

inline double EvalSinSum(const std::vector<double> &point) {
  double s = 0.0;
  for (double x : point) {
    s += std::sin(x);
  }
  return s;
}

inline double EvalSqrtSum(const std::vector<double> &point) {
  double s = 0.0;
  for (double x : point) {
    s += std::sqrt(x);
  }
  return s;
}

inline double EvalQuarticSum(const std::vector<double> &point) {
  double s = 0.0;
  for (double x : point) {
    s += x * x * x * x;
  }
  return s;
}

inline double EvaluateAt(FuncType func_type, const std::vector<double> &point) {
  switch (func_type) {
    case FuncType::kLinear:
      return EvalLinear(point);
    case FuncType::kSumCubes:
      return EvalSumCubes(point);
    case FuncType::kCosProduct:
      return EvalCosProduct(point);
    case FuncType::kExpNeg:
      return EvalExpNeg(point);
    case FuncType::kMixedPoly:
      return EvalMixedPoly(point);
    case FuncType::kSinSum:
      return EvalSinSum(point);
    case FuncType::kSqrtSum:
      return EvalSqrtSum(point);
    case FuncType::kQuarticSum:
      return EvalQuarticSum(point);
    default:
      return 0.0;
  }
}

}  // namespace detail

class SabirovSMonteCarloTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit SabirovSMonteCarloTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  std::vector<double> lower_;
  std::vector<double> upper_;
  int num_samples_{};
  FuncType func_type_{};
};

}  // namespace sabirov_s_monte_carlo
