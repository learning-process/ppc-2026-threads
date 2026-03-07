#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace batushin_i_incr_contrast_with_lhs {

using InType = std::vector<unsigned char>;
using OutType = std::vector<unsigned char>;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace batushin_i_incr_contrast_with_lhs
