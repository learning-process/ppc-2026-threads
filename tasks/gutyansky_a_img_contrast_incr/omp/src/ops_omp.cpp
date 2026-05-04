#include "gutyansky_a_img_contrast_incr/omp/include/ops_omp.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include "gutyansky_a_img_contrast_incr/common/include/common.hpp"

namespace gutyansky_a_img_contrast_incr {

GutyanskyAImgContrastIncrOMP::GutyanskyAImgContrastIncrOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool GutyanskyAImgContrastIncrOMP::ValidationImpl() {
  return !GetInput().empty();
}

bool GutyanskyAImgContrastIncrOMP::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return true;
}

bool GutyanskyAImgContrastIncrOMP::RunImpl() {
  const auto &input = GetInput();
  auto &output = GetOutput();

  const size_t sz = input.size();
  uint8_t lower_bound = std::numeric_limits<uint8_t>::max();
  uint8_t upper_bound = std::numeric_limits<uint8_t>::min();

#pragma omp parallel for default(none) shared(input, sz) reduction(min : lower_bound) reduction(max : upper_bound)
  for (size_t i = 0; i < sz; i++) {
    uint8_t val = input[i];
    lower_bound = std::min(lower_bound, val);
    upper_bound = std::max(upper_bound, val);
  }

  uint8_t delta = upper_bound - lower_bound;

  if (delta == 0) {
#pragma omp parallel for default(none) shared(input, output, sz)
    for (size_t i = 0; i < sz; i++) {
      output[i] = input[i];
    }
  } else {
#pragma omp parallel for default(none) shared(input, output, sz, lower_bound, delta)
    for (size_t i = 0; i < sz; i++) {
      auto old_value = static_cast<uint16_t>(input[i]);
      uint16_t new_value = (std::numeric_limits<uint8_t>::max() * (old_value - lower_bound)) / delta;

      output[i] = static_cast<uint8_t>(new_value);
    }
  }

  return true;
}

bool GutyanskyAImgContrastIncrOMP::PostProcessingImpl() {
  return true;
}

}  // namespace gutyansky_a_img_contrast_incr
