#include "pylaeva_s_inc_contrast_img_by_lsh/omp/include/ops_omp.hpp"

#include <algorithm>
#include <cmath>  // для std::round
#include <cstddef>
#include <cstdint>
#include <vector>

#include "pylaeva_s_inc_contrast_img_by_lsh/common/include/common.hpp"

namespace pylaeva_s_inc_contrast_img_by_lsh {

PylaevaSIncContrastImgByLshOMP::PylaevaSIncContrastImgByLshOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool PylaevaSIncContrastImgByLshOMP::ValidationImpl() {
  return !(GetInput().empty());
}

bool PylaevaSIncContrastImgByLshOMP::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return true;
}

bool PylaevaSIncContrastImgByLshOMP::RunImpl() {
  const auto &input = GetInput();
  auto &output = GetOutput();

  uint8_t min_pixel = input[0];
  uint8_t max_pixel = input[0];

  const size_t input_size = input.size();

#pragma omp parallel for default(none) shared(input, input_size) reduction(min : min_pixel) reduction(max : max_pixel)
  for (size_t i = 0; i < input_size; ++i) {
    min_pixel = std::min(min_pixel, input[i]);
    max_pixel = std::max(max_pixel, input[i]);
  }

  if (min_pixel == max_pixel) {
    output = input;
    return true;
  }

  float scale = 255.0F / static_cast<float>(max_pixel - min_pixel);

#pragma omp parallel for default(none) shared(input, output, min_pixel, input_size, scale)
  for (size_t i = 0; i < input.size(); ++i) {
    output[i] = static_cast<uint8_t>(std::round(static_cast<float>(input[i] - min_pixel) * scale));
  }

  return true;
}

bool PylaevaSIncContrastImgByLshOMP::PostProcessingImpl() {
  return true;
}

}  // namespace pylaeva_s_inc_contrast_img_by_lsh
