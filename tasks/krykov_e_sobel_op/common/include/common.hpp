#pragma once

#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace krykov_e_sobel_op {

struct Pixel {
  uint8_t r;
  uint8_t g;
  uint8_t b;

  bool operator==(const Pixel &other) const {
    return r == other.r && g == other.g && b == other.b;
  }
};

struct Image {
  int width = 0;
  int height = 0;
  std::vector<Pixel> data;
};

using InType = Image;
using OutType = std::vector<int>;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace krykov_e_sobel_op
