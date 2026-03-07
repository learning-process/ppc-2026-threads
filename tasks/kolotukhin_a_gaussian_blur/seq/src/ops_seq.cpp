#include "kolotukhin_a_gaussian_blur/seq/include/ops_seq.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "kolotukhin_a_gaussian_blur/common/include/common.hpp"

namespace kolotukhin_a_gaussian_blur {

KolotukhinAGaussinBlureSEQ::KolotukhinAGaussinBlureSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool KolotukhinAGaussinBlureSEQ::ValidationImpl() {
  const auto &pixel_data = get<0>(GetInput());
  const auto img_width = get<1>(GetInput());
  const auto img_height = get<2>(GetInput());

  return static_cast<std::size_t>(img_height) * static_cast<std::size_t>(img_width) == pixel_data.size();
}

bool KolotukhinAGaussinBlureSEQ::PreProcessingImpl() {
  const auto img_width = get<1>(GetInput());
  const auto img_height = get<2>(GetInput());

  GetOutput().assign(static_cast<std::size_t>(img_height) * static_cast<std::size_t>(img_width), 0);
  return true;
}

bool KolotukhinAGaussinBlureSEQ::RunImpl() {
  const auto &pixel_data = get<0>(GetInput());
  const auto img_width = get<1>(GetInput());
  const auto img_height = get<2>(GetInput());

  const static std::array<std::array<int, 3>, 3> kKernel = {{{{1, 2, 1}}, {{2, 4, 2}}, {{1, 2, 1}}}};
  const static int kSum = 16;

  auto &output = GetOutput();

  for (int row = 0; row < img_height; row++) {
    for (int col = 0; col < img_width; col++) {
      int acc = 0;
      for (int dy = -1; dy <= 1; dy++) {
        for (int dx = -1; dx <= 1; dx++) {
          std::uint8_t pixel = GetPixel(pixel_data, img_width, img_height, col + dx, row + dy);
          acc += kKernel.at(1 + dy).at(1 + dx) * static_cast<int>(pixel);
        }
      }
      output[(static_cast<std::size_t>(row) * static_cast<std::size_t>(img_width)) + static_cast<std::size_t>(col)] =
          static_cast<std::uint8_t>(acc / kSum);
    }
  }
  return true;
}

std::uint8_t KolotukhinAGaussinBlureSEQ::GetPixel(const std::vector<std::uint8_t> &pixel_data, int img_width,
                                                  int img_height, int pos_x, int pos_y) {
  std::size_t x = static_cast<std::size_t>(std::max(0, std::min(pos_x, img_width - 1)));
  std::size_t y = static_cast<std::size_t>(std::max(0, std::min(pos_y, img_height - 1)));
  return pixel_data[(y * static_cast<std::size_t>(img_width)) + x];
}

bool KolotukhinAGaussinBlureSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace kolotukhin_a_gaussian_blur
