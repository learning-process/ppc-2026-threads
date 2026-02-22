#include "kolotukhin_a_gaussian_blur/seq/include/ops_seq.hpp"

#include <numeric>
#include <vector>

#include "kolotukhin_a_gaussian_blur/common/include/common.hpp"
#include "util/include/util.hpp"

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

  return img_height * img_width == pixel_data.size();
}

bool KolotukhinAGaussinBlureSEQ::PreProcessingImpl() {
  const auto img_width = get<1>(GetInput());
  const auto img_height = get<2>(GetInput());

  GetOutput().assign(img_height * img_width, 0);
  return true;
}

bool KolotukhinAGaussinBlureSEQ::RunImpl() {
  const auto &pixel_data = get<0>(GetInput());
  const auto img_width = get<1>(GetInput());
  const auto img_height = get<2>(GetInput());

  const int kernel[3][3] = {{1, 2, 1}, {2, 4, 2}, {1, 2, 1}};
  const int kSum = 16;

  auto &output = GetOutput();

  for (int y = 0; y < img_height; y++) {
    for (int x = 0; x < img_width; x++) {
      int acc = 0;
      for (int dy = -1; dy <= 1; dy++) {
        for (int dx = -1; dx <= 1; dx++) {
          int pixel = Clamp(pixel_data, img_width, img_height, x + dx, y + dy);
          acc += kernel[1 + dy][1 + dx] * pixel;
        }
      }
      output[(static_cast<size_t>(y) * img_width) + static_cast<size_t>(x)] = static_cast<uint8_t>(acc / kSum);
    }
  }
  return true;
}

char8_t KolotukhinAGaussinBlureSEQ::Clamp(const std::vector<uint8_t> &pixel_data, int img_width, int img_height,
                                          int pos_x, int pos_y) const {
  std::size_t x = std::max(0, std::min(pos_x, img_width - 1));
  std::size_t y = std::max(0, std::min(pos_y, img_height - 1));
  return pixel_data[(y * static_cast<size_t>(img_width)) + x];
}

bool KolotukhinAGaussinBlureSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace kolotukhin_a_gaussian_blur
