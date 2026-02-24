#pragma once

#include <string>
#include <tuple>

#include "task/include/task.hpp"

namespace kamaletdinov_r_bitwise_int_seq {

using InType = int;
using OutType = int;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace kamaletdinov_r_bitwise_int_seq
