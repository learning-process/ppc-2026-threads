#pragma once

#include <vector>

#include "task/include/task.hpp"

namespace klimenko_v_lsh_contrast_incr_seq {

using InType = std::vector<int>;
using OutType = std::vector<int>;
using TestType = int;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace klimenko_v_lsh_contrast_incr_seq
