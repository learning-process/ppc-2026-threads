#pragma once

#include <vector>
#include "task/include/task.hpp"

namespace nazyrov_a_a_striped_multiplication {

using InType = std::pair<std::vector<double>, std::vector<double>>;  // A и B матрицы
using OutType = std::vector<double>;  // Результат C
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace nazyrov_a_a_striped_multiplication