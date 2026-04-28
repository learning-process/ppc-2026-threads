#include "batkov_f_contrast_enh_lin_hist_stretch/tbb/include/ops_tbb.hpp"

#include <oneapi/tbb.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

#include "batkov_f_contrast_enh_lin_hist_stretch/common/include/common.hpp"

namespace batkov_f_contrast_enh_lin_hist_stretch {

namespace {

std::pair<uint8_t, uint8_t> FindMinMaxParallel(const std::vector<uint8_t> &input) {
  using MinMax = std::pair<uint8_t, uint8_t>;
  const MinMax identity{std::numeric_limits<uint8_t>::max(), std::numeric_limits<uint8_t>::min()};

  return oneapi::tbb::parallel_reduce(oneapi::tbb::blocked_range<size_t>(0, input.size()), identity,
                                      [&](const oneapi::tbb::blocked_range<size_t> &r, MinMax mm) -> MinMax {
    for (size_t i = r.begin(); i != r.end(); ++i) {
      const uint8_t v = input[i];
      mm.first = std::min(mm.first, v);
      mm.second = std::max(mm.second, v);
    }
    return mm;
  }, [](const MinMax &a, const MinMax &b) -> MinMax {
    return {std::min(a.first, b.first), std::max(a.second, b.second)};
  });
}

std::pair<uint8_t, uint8_t> FindMinMax(const std::vector<uint8_t> &input, size_t parallel_threshold) {
  if (input.size() < parallel_threshold) {
    const auto [min_it, max_it] = std::ranges::minmax_element(input.begin(), input.end());
    return {*min_it, *max_it};
  }
  return FindMinMaxParallel(input);
}

}  // namespace

BatkovFContrastEnhLinHistStretchTBB::BatkovFContrastEnhLinHistStretchTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool BatkovFContrastEnhLinHistStretchTBB::ValidationImpl() {
  return !GetInput().empty();
}

bool BatkovFContrastEnhLinHistStretchTBB::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return true;
}

bool BatkovFContrastEnhLinHistStretchTBB::RunImpl() {
  auto &input = GetInput();
  auto &output = GetOutput();

  constexpr size_t kParallelMinMaxThreshold = 100000;
  const auto [min_el, max_el] = FindMinMax(input, kParallelMinMaxThreshold);

  if (max_el == min_el) {
    std::ranges::copy(input.begin(), input.end(), output.begin());
    return true;
  }

  const double a = 255.0 / static_cast<double>(max_el - min_el);
  const double b = -a * static_cast<double>(min_el);

  std::array<uint8_t, 256> lut{};
  for (size_t i = 0; i < lut.size(); ++i) {
    const double new_pixel = (a * static_cast<double>(i)) + b;
    lut.at(i) = static_cast<uint8_t>(std::clamp(new_pixel, 0.0, 255.0));
  }

  constexpr size_t kGrain = 1 << 17;
  oneapi::tbb::static_partitioner part;

  oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<size_t>(0, input.size(), kGrain),
                            [&](const oneapi::tbb::blocked_range<size_t> &r) {
    const uint8_t *src = input.data() + r.begin();
    uint8_t *dst = output.data() + r.begin();

    for (size_t i = 0; i < r.size(); ++i) {
      dst[i] = lut.at(src[i]);
    }
  }, part);

  return true;
}

bool BatkovFContrastEnhLinHistStretchTBB::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace batkov_f_contrast_enh_lin_hist_stretch
