#pragma once

#include <string>
#include <tuple>

#include "task/include/task.hpp"

namespace karpich_i_bitwise_batcher_seq {

using InType = int;
using OutType = int;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace karpich_i_bitwise_batcher_seq
