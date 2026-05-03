#pragma once

#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace balchunayte_z_sobel {

struct Pixel {
  uint8_t r = 0;
  uint8_t g = 0;
  uint8_t b = 0;
};

struct Image {
  int width = 0;
  int height = 0;
  std::vector<Pixel> data;  // size = width*height
};

using InType = Image;
using OutType = std::vector<int>;  // Sobel magnitude per pixel (abs(gx)+abs(gy))
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace balchunayte_z_sobel
