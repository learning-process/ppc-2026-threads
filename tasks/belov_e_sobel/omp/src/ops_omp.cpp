#include "belov_e_sobel/omp/include/ops_omp.hpp"

#include <omp.h>
#include <cmath>
#include <algorithm>
#include <vector>

#include "belov_e_sobel/common/include/common.hpp"
#include "util/include/util.hpp"

namespace belov_e_sobel {

BelovESobelOMP::BelovESobelOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = in;
}

bool BelovESobelOMP::ValidationImpl() {
  return !std::get<0>(GetInput()).empty() && (std::get<1>(GetInput()) > 0) && (std::get<2>(GetInput()) > 0);
}

bool BelovESobelOMP::PreProcessingImpl() {
  return true;
}

bool BelovESobelOMP::RunImpl() {
  const std::vector<uint8_t> &input = std::get<0>(GetInput());
  std::vector<uint8_t> &output = std::get<0>(GetOutput());
  int width = std::get<1>(GetInput());
  int height = std::get<2>(GetInput());

  auto get_px = [&](int x, int y) -> float {
    x = std::clamp(x, 0, width - 1);
    y = std::clamp(y, 0, height - 1);
    return static_cast<float>(input[(y * width) + x]);
  };

#pragma omp parallel for shared(input, output, width, height, get_px) schedule(dynamic)
  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      float gx = (-1 * get_px(x - 1, y - 1)) + (1 * get_px(x + 1, y - 1)) + (-2 * get_px(x - 1, y)) +
                 (2 * get_px(x + 1, y)) + (-1 * get_px(x - 1, y + 1)) + (1 * get_px(x + 1, y + 1));

      float gy = (-1 * get_px(x - 1, y - 1)) - (2 * get_px(x, y - 1)) - (1 * get_px(x + 1, y - 1)) +
                 (1 * get_px(x - 1, y + 1)) + (2 * get_px(x, y + 1)) + (1 * get_px(x + 1, y + 1));

      float magnitude = std::sqrt(gx * gx + gy * gy);
      output[y * width + x] = static_cast<uint8_t>(std::min(255.0f, magnitude));
    }
  }

  return true;
}

bool BelovESobelOMP::PostProcessingImpl() {
  return !std::get<0>(GetOutput()).empty() && (std::get<1>(GetOutput()) > 0) && (std::get<2>(GetOutput()) > 0);
}

}  // namespace belov_e_sobel
