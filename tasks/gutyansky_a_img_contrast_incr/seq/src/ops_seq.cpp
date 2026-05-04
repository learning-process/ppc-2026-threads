#include "gutyansky_a_img_contrast_incr/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>

#include "gutyansky_a_img_contrast_incr/common/include/common.hpp"

namespace gutyansky_a_img_contrast_incr {

GutyanskyAImgContrastIncrSEQ::GutyanskyAImgContrastIncrSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool GutyanskyAImgContrastIncrSEQ::ValidationImpl() {
  return !GetInput().empty();
}

bool GutyanskyAImgContrastIncrSEQ::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return true;
}

bool GutyanskyAImgContrastIncrSEQ::RunImpl() {
  const auto bounds = std::minmax_element(GetInput().begin(), GetInput().end());
  uint8_t lower_bound = *bounds.first;
  uint8_t upper_bound = *bounds.second;
  uint8_t delta = upper_bound - lower_bound;

  if (delta == 0) {
    std::copy(GetInput().begin(), GetInput().end(), GetOutput().begin());
  } else {
    size_t sz = GetInput().size();
    for (size_t i = 0; i < sz; i++) {
      auto old_value = static_cast<uint16_t>(GetInput()[i]);
      uint16_t new_value = (std::numeric_limits<uint8_t>::max() * (old_value - lower_bound)) / delta;

      GetOutput()[i] = static_cast<uint8_t>(new_value);
    }
  }

  return true;
}

bool GutyanskyAImgContrastIncrSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace gutyansky_a_img_contrast_incr
