#pragma once

#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace dorogin_v_bin_img_conv_hull_omp {

struct Point {
  int x{};
  int y{};
};

struct BinaryImage {
  int width{};
  int height{};
  std::vector<std::uint8_t> data;
};

using ComponentHull = std::vector<Point>;
using InType = BinaryImage;
using OutType = std::vector<ComponentHull>;

using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace dorogin_v_bin_img_conv_hull_omp
