#pragma once

#include <string>
#include <tuple>

#include "task/include/task.hpp"

namespace barkalova_m_mult_matrix_ccs {

using InType = int;
using OutType = int;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace barkalova_m_mult_matrix_ccs