#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace safronov_m_multiplication_matrix_blocksscheme_cannon {

using InType = std::tuple<int, std::vector<std::vector<double>>, std::vector<std::vector<double>>>;
using OutType = std::vector<std::vector<double>>;
using TestType = std::tuple<std::string, int, std::vector<std::vector<double>>, std::vector<std::vector<double>>,
                            std::vector<std::vector<double>>>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace safronov_m_multiplication_matrix_blocksscheme_cannon
