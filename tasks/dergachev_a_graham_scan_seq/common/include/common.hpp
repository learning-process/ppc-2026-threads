#pragma once

#include <string>
#include <tuple>

#include "task/include/task.hpp"

namespace dergachev_a_graham_scan_seq {

using InType = int;
using OutType = int;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace dergachev_a_graham_scan_seq
