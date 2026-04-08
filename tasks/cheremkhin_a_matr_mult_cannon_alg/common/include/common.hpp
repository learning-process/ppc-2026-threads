#pragma once

#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace cheremkhin_a_matr_mult_cannon_alg {

using InType = std::tuple<std::size_t, std::vector<double>, std::vector<double>>;
using OutType = std::vector<double>;
using TestType = std::tuple<std::size_t, int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace cheremkhin_a_matr_mult_cannon_alg
