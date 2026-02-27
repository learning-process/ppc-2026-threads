#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace pylaeva_s_inc_contrast_img_by_lsh {

using InType = std::vector<uint8_t>;
using OutType = std::vector<uint8_t>;
using TestType = std::tuple<std::vector<uint8_t>, std::vector<uint8_t>, std::string>;
;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace pylaeva_s_inc_contrast_img_by_lsh
