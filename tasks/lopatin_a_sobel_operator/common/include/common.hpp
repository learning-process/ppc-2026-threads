#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "task/include/task.hpp"

namespace lopatin_a_sobel_operator {

using PixelType = std::uint8_t;

struct Image {
  std::size_t height = 0;
  std::size_t width = 0;
  int threshold = 0;
  std::vector<PixelType> pixels;
};

using InType = Image;
using OutType = std::vector<int>;
using TestType = std::string;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace lopatin_a_sobel_operator
