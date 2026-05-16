#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace timofeev_n_radix_batcher_sort_threads {

using InType = std::vector<int>;
using OutType = std::vector<int>;
//                          вход,             выход,            проверка,  имя
using TestType = std::tuple<std::vector<int>, std::vector<int>, int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace timofeev_n_radix_batcher_sort_threads
