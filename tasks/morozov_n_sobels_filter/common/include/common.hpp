#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "task/include/task.hpp"

namespace morozov_n_sobels_filter {

struct Image {
  size_t height = 0;
  size_t width = 0;
  std::vector<uint8_t> pixels;
};

using InType = Image;
using OutType = Image;
using TestType = std::string;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace morozov_n_sobels_filter
