#pragma once

#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace fatehov_k_gaussian {

struct Image {
  uint32_t width = 0;
  uint32_t height = 0;
  uint32_t channels = 0;
  std::vector<uint8_t> data;

  Image() = default;
  Image(uint32_t w, uint32_t h, uint32_t ch) : width(w), height(h), channels(ch) {
    data.resize(static_cast<size_t>(w) * h * ch, 0);
  }
};

struct InputData {
  Image image;
  float sigma = 1.0f;
};

using InType = InputData;
using OutType = Image;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace fatehov_k_gaussian
