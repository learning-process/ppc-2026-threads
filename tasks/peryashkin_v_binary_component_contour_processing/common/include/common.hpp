// tasks/peryashkin_v_binary_component_contour_processing/common/include/common.hpp
#pragma once

#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace peryashkin_v_binary_component_contour_processing {

struct Point {
  int x{};
  int y{};
};

struct BinaryImage {
  int width{};
  int height{};
  std::vector<std::uint8_t> data;
};

using InType = BinaryImage;
using OutType = std::vector<std::vector<Point>>;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace peryashkin_v_binary_component_contour_processing
