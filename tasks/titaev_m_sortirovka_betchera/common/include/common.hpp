#pragma once
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace titaev_m_sortirovka_betchera {

using InType = std::vector<double>;
using OutType = std::vector<double>;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace titaev_m_sortirovka_betchera
