#include "pylaeva_s_inc_contrast_img_by_lsh/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>  // для std::round
#include <cstddef>
#include <cstdint>
#include <vector>

#include "pylaeva_s_inc_contrast_img_by_lsh/common/include/common.hpp"

namespace pylaeva_s_inc_contrast_img_by_lsh {

PylaevaSIncContrastImgByLshSEQ::PylaevaSIncContrastImgByLshSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool PylaevaSIncContrastImgByLshSEQ::ValidationImpl() {
  return !(GetInput().empty());
}

bool PylaevaSIncContrastImgByLshSEQ::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return true;
}

bool PylaevaSIncContrastImgByLshSEQ::RunImpl() {
  const auto &input = GetInput();
  auto &output = GetOutput();

  auto [min_pixel, max_pixel] = std::ranges::minmax(input);

  if (min_pixel == max_pixel) {
    output = input;
    return true;
  }

  float scale = 255.0F / static_cast<float>(max_pixel - min_pixel);

  for (size_t i = 0; i < input.size(); ++i) {
    output[i] = static_cast<uint8_t>(std::round(static_cast<float>(input[i] - min_pixel) * scale));
  }

  return true;
}

bool PylaevaSIncContrastImgByLshSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace pylaeva_s_inc_contrast_img_by_lsh
