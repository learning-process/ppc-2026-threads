#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace smyshlaev_a_sle_cg_seq {

using CGVector = std::vector<double>;
using CGMatrix = std::vector<CGVector>;

struct InputType {
  CGMatrix A;
  CGVector b;
};

using InType = InputType;
using OutType = CGVector;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace smyshlaev_a_sle_cg_seq
