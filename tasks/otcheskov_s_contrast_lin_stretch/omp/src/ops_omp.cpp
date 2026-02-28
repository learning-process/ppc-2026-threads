#include "otcheskov_s_contrast_lin_stretch/omp/include/ops_omp.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "otcheskov_s_contrast_lin_stretch/common/include/common.hpp"

namespace otcheskov_s_contrast_lin_stretch {

OtcheskovSContrastLinStretchOMP::OtcheskovSContrastLinStretchOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool OtcheskovSContrastLinStretchOMP::ValidationImpl() {
  return !GetInput().empty();
}

bool OtcheskovSContrastLinStretchOMP::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return GetOutput().size() == GetInput().size();
}

bool OtcheskovSContrastLinStretchOMP::RunImpl() {
  if (GetInput().empty()) {
    return false;
  }

  const InType &input = GetInput();
  OutType &output = GetOutput();

  uint8_t min_val = 255, max_val = 0;

#pragma omp parallel for reduction(min : min_val) reduction(max : max_val)
  for (size_t i = 0; i < input.size(); ++i) {
    uint8_t pixel = input[i];
    if (pixel < min_val) {
      min_val = pixel;
    }
    if (pixel > max_val) {
      max_val = pixel;
    }
  }

  const size_t min_copy_size = 1000000;
  if (min_val == max_val) {
#pragma omp parallel for if (input.size() > min_copy_size)
    for (size_t i = 0; i < input.size(); ++i) {
      output[i] = input[i];
    }

    return true;
  }

  const int min_i = static_cast<int>(min_val);
  const int range = static_cast<int>(max_val) - min_i;

#pragma omp parallel for
  for (size_t i = 0; i < input.size(); ++i) {
    int pixel = static_cast<int>(input[i]);
    int value = (pixel - min_i) * 255 / range;
    output[i] = static_cast<uint8_t>(std::clamp(value, 0, 255));
  }

  return true;
}

bool OtcheskovSContrastLinStretchOMP::PostProcessingImpl() {
  return true;
}

}  // namespace otcheskov_s_contrast_lin_stretch
