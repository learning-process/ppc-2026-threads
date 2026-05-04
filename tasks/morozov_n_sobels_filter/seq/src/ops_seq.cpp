#include "morozov_n_sobels_filter/seq/include/ops_seq.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>

#include "morozov_n_sobels_filter/common/include/common.hpp"

namespace morozov_n_sobels_filter {

MorozovNSobelsFilterSEQ::MorozovNSobelsFilterSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;

  result_image_.height = in.height;
  result_image_.width = in.width;
  result_image_.pixels.resize(result_image_.height * result_image_.width, 0);
}

bool MorozovNSobelsFilterSEQ::ValidationImpl() {
  const Image &input = GetInput();
  return (input.height == result_image_.height) && (input.width == result_image_.width) &&
         (input.pixels.size() == result_image_.pixels.size());
}

bool MorozovNSobelsFilterSEQ::PreProcessingImpl() {
  return true;
}

bool MorozovNSobelsFilterSEQ::RunImpl() {
  const Image &input = GetInput();
  Filter(input);
  GetOutput() = result_image_;
  return true;
}

bool MorozovNSobelsFilterSEQ::PostProcessingImpl() {
  return true;
}

void MorozovNSobelsFilterSEQ::Filter(const Image &img) {
  for (size_t id_y = 1; id_y < img.height - 1; id_y++) {
    for (size_t id_x = 1; id_x < img.width - 1; id_x++) {
      size_t pixel_id = (id_y * img.width) + id_x;
      result_image_.pixels[pixel_id] = CalculateNewPixelColor(img, id_x, id_y);
    }
  }
}

uint8_t MorozovNSobelsFilterSEQ::CalculateNewPixelColor(const Image &img, size_t x, size_t y) {
  constexpr int kRadX = 1;
  constexpr int kRadY = 1;
  constexpr size_t kZero = 0;

  int grad_x = 0;
  int grad_y = 0;

  for (int row_offset = -kRadY; row_offset <= kRadY; row_offset++) {
    for (int col_offset = -kRadX; col_offset <= kRadX; col_offset++) {
      size_t id_x = std::clamp(x + col_offset, kZero, img.width - 1);
      size_t id_y = std::clamp(y + row_offset, kZero, img.height - 1);
      size_t pixel_id = (id_y * img.width) + id_x;

      grad_x += img.pixels[pixel_id] * kKernelX.at(row_offset + kRadY).at(col_offset + kRadX);
      grad_y += img.pixels[pixel_id] * kKernelY.at(row_offset + kRadY).at(col_offset + kRadX);
    }
  }

  int gradient = static_cast<int>(std::sqrt((grad_x * grad_x) + (grad_y * grad_y)));
  gradient = std::clamp(gradient, 0, 255);

  return static_cast<uint8_t>(gradient);
}
}  // namespace morozov_n_sobels_filter
