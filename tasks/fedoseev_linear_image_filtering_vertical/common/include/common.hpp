#pragma once

#include <string>
#include <tuple>

#include "task/include/task.hpp"

namespace fedoseev_linear_image_filtering_vertical {

struct Image {
    int width;
    int height;
    std::vector<int> data;
};

using InType = Image;
using OutType = Image;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace fedoseev_linear_image_filtering_vertical