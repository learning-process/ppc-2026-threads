#include "zaharov_g_linear_contrast_stretch/tbb/include/ops_tbb.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>

#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/parallel_for.h"
#include "oneapi/tbb/parallel_reduce.h"
#include "zaharov_g_linear_contrast_stretch/common/include/common.hpp"

namespace zaharov_g_linear_contrast_stretch {

namespace {

struct MinMax {
  int min;
  int max;
};

MinMax MergeMinMax(const MinMax &left, const MinMax &right) {
  return MinMax{.min = std::min(left.min, right.min), .max = std::max(left.max, right.max)};
}

MinMax FindMinMax(const InType &input) {
  const auto range = oneapi::tbb::blocked_range<std::size_t>(0, input.size());
  const auto identity = MinMax{.min = std::numeric_limits<uint8_t>::max(), .max = std::numeric_limits<uint8_t>::min()};

  const auto find_local_minmax = [&input](const oneapi::tbb::blocked_range<std::size_t> &range, MinMax current) {
    for (std::size_t i = range.begin(); i != range.end(); ++i) {
      const int value = static_cast<int>(input[i]);
      current.min = std::min(current.min, value);
      current.max = std::max(current.max, value);
    }
    return current;
  };

  const auto merge_minmax = [](const MinMax &left, const MinMax &right) { return MergeMinMax(left, right); };

  return oneapi::tbb::parallel_reduce(range, identity, find_local_minmax, merge_minmax);
}

void StretchImage(const InType &input, OutType &output, const MinMax &minmax) {
  if (minmax.max > minmax.min) {
    const int denom = minmax.max - minmax.min;
    oneapi::tbb::parallel_for(
        oneapi::tbb::blocked_range<std::size_t>(0, input.size()),
        [&input, &output, min_el = minmax.min, denom](const oneapi::tbb::blocked_range<std::size_t> &range) {
      for (std::size_t i = range.begin(); i != range.end(); ++i) {
        const int value = (static_cast<int>(input[i]) - min_el) * std::numeric_limits<uint8_t>::max() / denom;
        output[i] = static_cast<uint8_t>(std::clamp(value, 0, 255));
      }
    });
  } else {
    oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<std::size_t>(0, input.size()),
                              [&input, &output](const oneapi::tbb::blocked_range<std::size_t> &range) {
      std::copy(input.begin() + static_cast<std::ptrdiff_t>(range.begin()),
                input.begin() + static_cast<std::ptrdiff_t>(range.end()),
                output.begin() + static_cast<std::ptrdiff_t>(range.begin()));
    });
  }
}

}  // namespace

ZaharovGLinContrStrTBB::ZaharovGLinContrStrTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool ZaharovGLinContrStrTBB::ValidationImpl() {
  return !GetInput().empty();
}

bool ZaharovGLinContrStrTBB::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return true;
}

bool ZaharovGLinContrStrTBB::RunImpl() {
  const InType &input = GetInput();
  OutType &output = GetOutput();

  const MinMax minmax = FindMinMax(input);
  StretchImage(input, output, minmax);

  return true;
}

bool ZaharovGLinContrStrTBB::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace zaharov_g_linear_contrast_stretch
