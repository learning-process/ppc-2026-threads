#pragma once

#include <functional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "task/include/task.hpp"

namespace kutergin_v_multidimensional_integration_rect_method {

struct IntegrationData {
  std::vector<std::pair<double, double>> limits;            // вектор пар {a, b} - границ каждой из n размерностей
  std::vector<int> n_steps;                                 // количество шагов для каждой из n размерностей
  std::function<double(const std::vector<double> &)> func;  // сама функция
};

using InType = IntegrationData;
using OutType = double;
using BaseTask = ppc::task::Task<InType, OutType>;
using TestType = std::tuple<InType, OutType, std::string>;

}  // namespace kutergin_v_multidimensional_integration_rect_method
