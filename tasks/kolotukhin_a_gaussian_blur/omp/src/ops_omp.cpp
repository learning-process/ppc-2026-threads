#include "kolotukhin_a_gaussian_blur/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <tuple>
#include <vector>

#include "kolotukhin_a_gaussian_blur/common/include/common.hpp"

namespace kolotukhin_a_gaussian_blur {

KolotukhinAGaussinBlureOMP::KolotukhinAGaussinBlureOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool KolotukhinAGaussinBlureOMP::ValidationImpl() {
  const auto &pixel_data = std::get<0>(GetInput());
  const auto img_width = std::get<1>(GetInput());
  const auto img_height = std::get<2>(GetInput());

  return static_cast<std::size_t>(img_height) * static_cast<std::size_t>(img_width) == pixel_data.size();
}

bool KolotukhinAGaussinBlureOMP::PreProcessingImpl() {
  const auto img_width = std::get<1>(GetInput());
  const auto img_height = std::get<2>(GetInput());

  GetOutput().assign(static_cast<std::size_t>(img_height) * static_cast<std::size_t>(img_width), 0);
  return true;
}

bool KolotukhinAGaussinBlureOMP::RunImpl() {
  const auto &pixel_data = std::get<0>(GetInput());
  const auto img_width = std::get<1>(GetInput());
  const auto img_height = std::get<2>(GetInput());

  const static std::array<int, 9> kKernel = {1, 2, 1, 2, 4, 2, 1, 2, 1};
  const static int kSum = 16;

  auto &output = GetOutput();

#pragma omp parallel for collapse(2) schedule(static) default(none) \
    shared(pixel_data, img_width, img_height, output, kKernel, kSum)
  for (int row = 0; row < img_height; row++) {
    for (int col = 0; col < img_width; col++) {
      int acc = 0;
      for (int dy = -1; dy <= 1; dy++) {
        for (int dx = -1; dx <= 1; dx++) {
          std::uint8_t pixel = GetPixel(pixel_data, img_width, img_height, col + dx, row + dy);
          int idx = ((dy + 1) * 3) + (dx + 1);
          acc += kKernel.at(idx) * static_cast<int>(pixel);
        }
      }
      output[(static_cast<std::size_t>(row) * static_cast<std::size_t>(img_width)) + static_cast<std::size_t>(col)] =
          static_cast<std::uint8_t>(acc / kSum);
    }
  }
  return true;
}

std::uint8_t KolotukhinAGaussinBlureOMP::GetPixel(const std::vector<std::uint8_t> &pixel_data, int img_width,
                                                  int img_height, int pos_x, int pos_y) {
  std::size_t x = static_cast<std::size_t>(std::max(0, std::min(pos_x, img_width - 1)));
  std::size_t y = static_cast<std::size_t>(std::max(0, std::min(pos_y, img_height - 1)));
  return pixel_data[(y * static_cast<std::size_t>(img_width)) + x];
}

bool KolotukhinAGaussinBlureOMP::PostProcessingImpl() {
  return true;
}

}  // namespace kolotukhin_a_gaussian_blur
