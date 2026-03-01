#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace chaschin_v_linear_image_filtration_seq {

using InType = std::tuple<std::vector<float>,int,int>;
using OutType = std::vector<float>;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace chaschin_v_linear_image_filtration_seq
