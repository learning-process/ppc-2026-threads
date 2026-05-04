#pragma once

#include <functional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "task/include/task.hpp"

namespace chernykh_s_trapezoidal_integration {

struct IntegrationInType {
  std::vector<std::pair<double, double>> limits;
  std::vector<int> steps;
  std::function<double(const std::vector<double> &)> func;

  IntegrationInType() = default;

  IntegrationInType(std::vector<std::pair<double, double>> l, std::vector<int> s,
                    std::function<double(const std::vector<double> &)> f)
      : limits(std::move(l)), steps(std::move(s)), func(std::move(f)) {}
};

using InType = IntegrationInType;
using OutType = double;
using TestType = std::tuple<InType, OutType, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace chernykh_s_trapezoidal_integration
