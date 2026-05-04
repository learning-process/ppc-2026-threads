#pragma once

#include <vector>

#include "task/include/task.hpp"

namespace lukin_i_ench_contr_lin_hist {

using InType = std::vector<unsigned char>;
using OutType = std::vector<unsigned char>;
using TestType = int;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace lukin_i_ench_contr_lin_hist
