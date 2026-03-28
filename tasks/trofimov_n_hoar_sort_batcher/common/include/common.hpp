#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace trofimov_n_hoar_sort_batcher {

using InType = std::vector<int>;
using OutType = std::vector<int>;
using TestType = std::tuple<std::vector<int>, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace trofimov_n_hoar_sort_batcher
