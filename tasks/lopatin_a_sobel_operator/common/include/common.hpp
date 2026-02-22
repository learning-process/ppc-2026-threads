#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace lopatin_a_sobel_operator {

using PixelType = std::uint8_t;

struct Image {
  std::size_t height;
  std::size_t width;
  int threshold;
  std::vector<PixelType> pixels;
};

using InType = Image;
using OutType = std::vector<PixelType>;
using TestType = std::string;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace lopatin_a_sobel_operator
