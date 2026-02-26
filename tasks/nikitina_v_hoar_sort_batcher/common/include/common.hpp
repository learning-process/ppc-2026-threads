#pragma once

#include <string>
#include <vector>
#include <tuple>


#include "task/include/task.hpp"

namespace nikitina_v_hoar_sort_batcher {

using InType = std::vector<int>;
using OutType = std::vector<int>;
using TestType = std::tuple<int, std::string>;

using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace nikitina_v_hoar_sort_batcher
