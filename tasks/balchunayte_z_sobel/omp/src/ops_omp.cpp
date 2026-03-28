#include "balchunayte_z_sobel/omp/include/ops_omp.hpp"

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

BalchunayteZSobelOpOMP::BalchunayteZSobelOpOMP(const InType &input_image) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = input_image;
  GetOutput().clear();
}

bool BalchunayteZSobelOpOMP::ValidationImpl() {
  const auto &input_image = GetInput();

  if (input_image.width <= 0 || input_image.height <= 0) {
    return false;
  }

  const auto expected_size = static_cast<size_t>(input_image.width) * static_cast<size_t>(input_image.height);

  if (input_image.data.size() != expected_size) {
    return false;
  }

  return GetOutput().empty();
}

bool BalchunayteZSobelOpOMP::PreProcessingImpl() {
  const auto &input_image = GetInput();
  GetOutput().assign(static_cast<size_t>(input_image.width) * static_cast<size_t>(input_image.height), 0);
  return true;
}

bool BalchunayteZSobelOpOMP::RunImpl() {
  const auto &input_image = GetInput();
  auto &output_data = GetOutput();

  const int image_width = input_image.width;
  const int image_height = input_image.height;
  const auto image_width_size = static_cast<size_t>(image_width);

  if (image_width < 3 || image_height < 3) {
    return true;
  }

#pragma omp parallel for default(none) shared(input_image, output_data, image_width, image_height, image_width_size) \
    schedule(static)
  for (int row_index = 1; row_index < image_height - 1; ++row_index) {
    for (int col_index = 1; col_index < image_width - 1; ++col_index) {
      const size_t index_top_left =
          (static_cast<size_t>(row_index - 1) * image_width_size) + static_cast<size_t>(col_index - 1);
      const size_t index_top_middle =
          (static_cast<size_t>(row_index - 1) * image_width_size) + static_cast<size_t>(col_index);
      const size_t index_top_right =
          (static_cast<size_t>(row_index - 1) * image_width_size) + static_cast<size_t>(col_index + 1);

      const size_t index_middle_left =
          (static_cast<size_t>(row_index) * image_width_size) + static_cast<size_t>(col_index - 1);
      const size_t index_middle_right =
          (static_cast<size_t>(row_index) * image_width_size) + static_cast<size_t>(col_index + 1);

      const size_t index_bottom_left =
          (static_cast<size_t>(row_index + 1) * image_width_size) + static_cast<size_t>(col_index - 1);
      const size_t index_bottom_middle =
          (static_cast<size_t>(row_index + 1) * image_width_size) + static_cast<size_t>(col_index);
      const size_t index_bottom_right =
          (static_cast<size_t>(row_index + 1) * image_width_size) + static_cast<size_t>(col_index + 1);

      const int gray_top_left = ConvertPixelToGray(input_image.data[index_top_left]);
      const int gray_top_middle = ConvertPixelToGray(input_image.data[index_top_middle]);
      const int gray_top_right = ConvertPixelToGray(input_image.data[index_top_right]);

      const int gray_middle_left = ConvertPixelToGray(input_image.data[index_middle_left]);
      const int gray_middle_right = ConvertPixelToGray(input_image.data[index_middle_right]);

      const int gray_bottom_left = ConvertPixelToGray(input_image.data[index_bottom_left]);
      const int gray_bottom_middle = ConvertPixelToGray(input_image.data[index_bottom_middle]);
      const int gray_bottom_right = ConvertPixelToGray(input_image.data[index_bottom_right]);

      const int gradient_x = (-gray_top_left + gray_top_right) + (-2 * gray_middle_left + 2 * gray_middle_right) +
                             (-gray_bottom_left + gray_bottom_right);

      const int gradient_y = (gray_top_left + (2 * gray_top_middle) + gray_top_right) +
                             (-gray_bottom_left - (2 * gray_bottom_middle) - gray_bottom_right);

      const int magnitude = std::abs(gradient_x) + std::abs(gradient_y);

      const size_t output_index = (static_cast<size_t>(row_index) * image_width_size) + static_cast<size_t>(col_index);
      output_data[output_index] = magnitude;
    }
  }

  return true;
}

bool BalchunayteZSobelOpOMP::PostProcessingImpl() {
  return true;
}

}  // namespace balchunayte_z_sobel
