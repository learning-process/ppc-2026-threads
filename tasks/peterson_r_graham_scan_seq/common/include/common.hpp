#pragma once

#include <string>
#include <tuple>

#include "task/include/task.hpp"

namespace peterson_r_graham_scan_seq {

using InputValue = int;
using OutputValue = int;
using TestParameters = std::tuple<int, std::string>;
using TaskBase = ppc::task::Task<InputValue, OutputValue>;

}  // namespace peterson_r_graham_scan_seq
