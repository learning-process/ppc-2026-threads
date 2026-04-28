#include "gutyansky_a_img_contrast_incr/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>

#include "gutyansky_a_img_contrast_incr/common/include/common.hpp"
#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/parallel_for.h"

namespace gutyansky_a_img_contrast_incr {

GutyanskyAImgContrastIncrTBB::GutyanskyAImgContrastIncrTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool GutyanskyAImgContrastIncrTBB::ValidationImpl() {
  return !GetInput().empty();
}

bool GutyanskyAImgContrastIncrTBB::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return true;
}

bool GutyanskyAImgContrastIncrTBB::RunImpl() {
  const auto &input = GetInput();
  auto &output = GetOutput();

  const size_t sz = input.size();
  auto [lower_bound, upper_bound] = tbb::parallel_reduce(
      tbb::blocked_range<size_t>(0, sz), std::make_pair(static_cast<uint8_t>(255), static_cast<uint8_t>(0)),
      [&input](const auto &range, auto init) {
    auto [local_min, local_max] = init;
    for (size_t i = range.begin(); i != range.end(); ++i) {
      local_min = std::min(local_min, input[i]);
      local_max = std::max(local_max, input[i]);
    }
    return std::make_pair(local_min, local_max);
  }, [](auto a, auto b) { return std::make_pair(std::min(a.first, b.first), std::max(a.second, b.second)); });

  uint8_t delta = upper_bound - lower_bound;

  if (delta == 0) {
    tbb::parallel_for(tbb::blocked_range<size_t>(0, sz), [&](const auto &range) {
      for (auto idx = range.begin(); idx != range.end(); ++idx) {
        output[idx] = input[idx];
      }
    });
  } else {
    constexpr uint16_t kMaxUint8 = std::numeric_limits<uint8_t>::max();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, sz), [&](const auto &range) {
      for (auto idx = range.begin(); idx != range.end(); ++idx) {
        uint16_t old_value = input[idx];
        uint16_t new_value = (kMaxUint8 * (old_value - lower_bound)) / delta;
        output[idx] = static_cast<uint8_t>(new_value);
      }
    });
  }

  return true;
}

bool GutyanskyAImgContrastIncrTBB::PostProcessingImpl() {
  return true;
}

}  // namespace gutyansky_a_img_contrast_incr
