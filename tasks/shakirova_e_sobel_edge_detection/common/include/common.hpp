#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"
#include "img_container.hpp"

namespace shakirova_e_sobel_edge_detection {

using InType   = ImgContainer;
using OutType  = std::vector<int>;   
using TestType = std::tuple<int, std::string>;  
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace shakirova_e_sobel_edge_detection