#include "batkov_f_contrast_enh_lin_hist_stretch_seq/omp/include/ops_omp.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

#include "batkov_f_contrast_enh_lin_hist_stretch_seq/common/include/common.hpp"

namespace {

std::pair<uint8_t, uint8_t> FindMinMaxParallel(const std::vector<uint8_t> &input) {
  uint8_t min_el = std::numeric_limits<uint8_t>::max();
  uint8_t max_el = std::numeric_limits<uint8_t>::min();

#pragma omp parallel for default(none) shared(input) reduction(min : min_el) reduction(max : max_el)
  for (uint8_t p : input) {
    min_el = std::min(p, min_el);
    max_el = std::max(p, max_el);
  }

  return { min_el, max_el };
}

std::pair<uint8_t, uint8_t> FindMinMax(const std::vector<uint8_t> &input, size_t parallel_threshold) {
  if (input.size() < parallel_threshold) {
    auto [it_min, it_max] = std::ranges::minmax_element(input);
    return {*it_min, *it_max};
  }

  return FindMinMaxParallel(input);
}

}  // namespace

namespace batkov_f_contrast_enh_lin_hist_stretch_seq {

BatkovFContrastEnhLinHistStretchOMP::BatkovFContrastEnhLinHistStretchOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool BatkovFContrastEnhLinHistStretchOMP::ValidationImpl() {
  return !GetInput().empty();
}

bool BatkovFContrastEnhLinHistStretchOMP::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return true;
}

bool BatkovFContrastEnhLinHistStretchOMP::RunImpl() {
  const auto &input = GetInput();
  auto &output = GetOutput();

  auto [min_el, max_el] = FindMinMax(input, 100000);

  if (max_el > min_el) {
    double a = 255.0 / (max_el - min_el);
    double b = -a * min_el;

#pragma omp parallel for default(none) shared(input, output, a, b)
    for (size_t i = 0; i < input.size(); ++i) {
      double new_pixel = (a * static_cast<double>(input[i])) + b;
      output[i] = static_cast<uint8_t>(std::clamp(new_pixel, 0.0, 255.0));
    }
  } else {
    std::ranges::copy(input.begin(), input.end(), output.begin());
  }

  return true;
}

bool BatkovFContrastEnhLinHistStretchOMP::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace batkov_f_contrast_enh_lin_hist_stretch_seq
