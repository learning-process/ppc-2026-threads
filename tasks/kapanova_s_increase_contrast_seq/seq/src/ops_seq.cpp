#include "kapanova_s_increase_contrast_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace kapanova_s_increase_contrast_seq {

KapanovaSIncreaseContrastSEQ::KapanovaSIncreaseContrastSEQ(const InType &in) : min_pixel_(0), max_pixel_(0) {
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
    std::ranges::fill(output, min_pixel_);
    return true;
  }

  const float scale = 255.0F / static_cast<float>(max_pixel_ - min_pixel_);

  for (size_t i = 0; i < input.size(); ++i) {
    float normalized = static_cast<float>(input[i] - min_pixel_) * scale;
    output[i] = static_cast<uint8_t>(std::lround(normalized));
  }

  return true;
}

bool KapanovaSIncreaseContrastSEQ::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace kapanova_s_increase_contrast_seq
