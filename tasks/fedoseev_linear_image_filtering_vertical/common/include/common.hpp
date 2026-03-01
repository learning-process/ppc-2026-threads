#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace fedoseev_linear_image_filtering_vertical {

struct Image {
  int width = 0;
  int height = 0;
  std::vector<int> data;

  Image() = default;
  Image(int w, int h, const std::vector<int> &d) : width(w), height(h), data(d) {}
};

using InType = Image;
using OutType = Image;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace fedoseev_linear_image_filtering_vertical
