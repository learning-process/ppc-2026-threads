#pragma once

#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace makoveeva_matmul_double_seq {

using InType = std::tuple<size_t, std::vector<double>, std::vector<double>>;
using OutType = std::vector<double>;
using TestType = std::tuple<size_t, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace makoveeva_matmul_double_seq
