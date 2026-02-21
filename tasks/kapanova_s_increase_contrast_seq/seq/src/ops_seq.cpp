#include "kapanova_s_increase_contrast_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstdint>
#include <vector>

namespace kapanova_s_increase_contrast_seq {

KapanovaSIncreaseContrastSEQ::KapanovaSIncreaseContrastSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType();
}

bool KapanovaSIncreaseContrastSEQ::ValidationImpl() {
  return !GetInput().empty();
}

bool KapanovaSIncreaseContrastSEQ::PreProcessingImpl() {
  const auto &input = GetInput();
  min_pixel_ = *std::min_element(input.begin(), input.end());
  max_pixel_ = *std::max_element(input.begin(), input.end());
  return (max_pixel_ >= min_pixel_);
}

bool KapanovaSIncreaseContrastSEQ::RunImpl() {
  const auto &input = GetInput();
  auto &output = GetOutput();
  output.resize(input.size());

  if (max_pixel_ == min_pixel_) {
    std::fill(output.begin(), output.end(), min_pixel_);
    return true;
  }

  const float scale = 255.0f / static_cast<float>(max_pixel_ - min_pixel_);

  for (size_t i = 0; i < input.size(); ++i) {
    float normalized = static_cast<float>(input[i] - min_pixel_) * scale;
    output[i] = static_cast<uint8_t>(normalized + 0.5f);
  }

  return true;
}

bool KapanovaSIncreaseContrastSEQ::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace kapanova_s_increase_contrast_seq
