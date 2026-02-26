#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "img_container.hpp"
#include "task/include/task.hpp"

namespace shakirova_e_sobel_edge_detection {

using InType   = ImgContainer;
using OutType  = std::vector<int>;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace shakirova_e_sobel_edge_detection