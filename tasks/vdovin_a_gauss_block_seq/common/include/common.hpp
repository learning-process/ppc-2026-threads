#pragma once

#include <string>
#include <tuple>

#include "task/include/task.hpp"

namespace vdovin_a_gauss_block_seq {

using InType = int;
using OutType = int;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace vdovin_a_gauss_block_seq
