#include "otcheskov_s_contrast_lin_stretch/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "otcheskov_s_contrast_lin_stretch/common/include/common.hpp"

namespace otcheskov_s_contrast_lin_stretch {

OtcheskovSContrastLinStretchSEQ::OtcheskovSContrastLinStretchSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool OtcheskovSContrastLinStretchSEQ::ValidationImpl() {
  return !GetInput().empty();
}

bool OtcheskovSContrastLinStretchSEQ::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return GetOutput().size() == GetInput().size();
}

bool OtcheskovSContrastLinStretchSEQ::RunImpl() {
  if (GetInput().empty()) {
    return false;
  }

  const InType &input = GetInput();
  OutType &output = GetOutput();

  auto [min_it, max_it] = std::ranges::minmax_element(input);
  uint8_t min_val = *min_it;
  uint8_t max_val = *max_it;

  if (min_val == max_val) {
    output = input;
    return true;
  }

  const int min_i = static_cast<int>(min_val);
  const int max_i = static_cast<int>(max_val);
  const int range = max_i - min_i;
  for (size_t i = 0; i < input.size(); ++i) {
    int pixel = static_cast<int>(input[i]);
    int value = (pixel - min_i) * 255 / (range);

    value = std::clamp(value, 0, 255);

    output[i] = static_cast<uint8_t>(value);
  }

  return true;
}

bool OtcheskovSContrastLinStretchSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace otcheskov_s_contrast_lin_stretch
