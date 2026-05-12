#include "balchunayte_z_sobel/seq/include/ops_seq.hpp"

#include <cstddef>
#include <cstdlib>
#include <vector>

#include "balchunayte_z_sobel/common/include/common.hpp"

namespace balchunayte_z_sobel {

namespace {

int ConvertPixelToGray(const Pixel &pixel_value) {
  return (77 * static_cast<int>(pixel_value.r) + 150 * static_cast<int>(pixel_value.g) +
          29 * static_cast<int>(pixel_value.b)) >>
         8;
}

}  // namespace

BalchunayteZSobelOpSEQ::BalchunayteZSobelOpSEQ(const InType &input_image) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = input_image;
  GetOutput().clear();
}

bool BalchunayteZSobelOpSEQ::ValidationImpl() {
  const auto &input_image = GetInput();

  if (input_image.width <= 0 || input_image.height <= 0) {
    return false;
  }

  const size_t expected_size = static_cast<size_t>(input_image.width) * static_cast<size_t>(input_image.height);

  if (input_image.data.size() != expected_size) {
    return false;
  }

  return GetOutput().empty();
}

bool BalchunayteZSobelOpSEQ::PreProcessingImpl() {
  const auto &input_image = GetInput();
  GetOutput().assign(static_cast<size_t>(input_image.width) * static_cast<size_t>(input_image.height), 0);
  return true;
}

bool BalchunayteZSobelOpSEQ::RunImpl() {
  const auto &input_image = GetInput();

  const int image_width = input_image.width;
  const int image_height = input_image.height;

  if (image_width < 3 || image_height < 3) {
    return true;
  }

  for (int row_index = 1; row_index < image_height - 1; ++row_index) {
    for (int col_index = 1; col_index < image_width - 1; ++col_index) {
      const int index_top_left = ((row_index - 1) * image_width) + (col_index - 1);
      const int index_top_middle = ((row_index - 1) * image_width) + (col_index + 0);
      const int index_top_right = ((row_index - 1) * image_width) + (col_index + 1);

      const int index_middle_left = ((row_index + 0) * image_width) + (col_index - 1);
      const int index_middle_right = ((row_index + 0) * image_width) + (col_index + 1);

      const int index_bottom_left = ((row_index + 1) * image_width) + (col_index - 1);
      const int index_bottom_middle = ((row_index + 1) * image_width) + (col_index + 0);
      const int index_bottom_right = ((row_index + 1) * image_width) + (col_index + 1);

      const int gray_top_left = ConvertPixelToGray(input_image.data[static_cast<size_t>(index_top_left)]);
      const int gray_top_middle = ConvertPixelToGray(input_image.data[static_cast<size_t>(index_top_middle)]);
      const int gray_top_right = ConvertPixelToGray(input_image.data[static_cast<size_t>(index_top_right)]);

      const int gray_middle_left = ConvertPixelToGray(input_image.data[static_cast<size_t>(index_middle_left)]);
      const int gray_middle_right = ConvertPixelToGray(input_image.data[static_cast<size_t>(index_middle_right)]);

      const int gray_bottom_left = ConvertPixelToGray(input_image.data[static_cast<size_t>(index_bottom_left)]);
      const int gray_bottom_middle = ConvertPixelToGray(input_image.data[static_cast<size_t>(index_bottom_middle)]);
      const int gray_bottom_right = ConvertPixelToGray(input_image.data[static_cast<size_t>(index_bottom_right)]);

      const int gradient_x = (-gray_top_left + gray_top_right) + (-2 * gray_middle_left + 2 * gray_middle_right) +
                             (-gray_bottom_left + gray_bottom_right);

      const int gradient_y = (gray_top_left + 2 * gray_top_middle + gray_top_right) +
                             (-gray_bottom_left - 2 * gray_bottom_middle - gray_bottom_right);

      const int magnitude = std::abs(gradient_x) + std::abs(gradient_y);

      const size_t output_index =
          (static_cast<size_t>(row_index) * static_cast<size_t>(image_width)) + static_cast<size_t>(col_index);
      GetOutput()[output_index] = magnitude;
    }
  }

  return true;
}

bool BalchunayteZSobelOpSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace balchunayte_z_sobel
