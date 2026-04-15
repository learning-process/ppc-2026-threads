#include "batkov_f_contrast_enh_lin_hist_stretch/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "batkov_f_contrast_enh_lin_hist_stretch/common/include/common.hpp"

namespace batkov_f_contrast_enh_lin_hist_stretch {

BatkovFContrastEnhLinHistStretchSEQ::BatkovFContrastEnhLinHistStretchSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool BatkovFContrastEnhLinHistStretchSEQ::ValidationImpl() {
  return !GetInput().empty();
}

bool BatkovFContrastEnhLinHistStretchSEQ::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return true;
}

bool BatkovFContrastEnhLinHistStretchSEQ::RunImpl() {
  auto &input = GetInput();
  auto &output = GetOutput();

  auto [min_it, max_it] = std::ranges::minmax_element(GetInput());
  uint8_t min_el = *min_it;
  uint8_t max_el = *max_it;

  if (max_el > min_el) {
    double a = 255.0 / (max_el - min_el);
    double b = -a * min_el;

    for (size_t i = 0; i < input.size(); ++i) {
      double new_pixel = (a * static_cast<double>(input[i])) + b;
      output[i] = static_cast<uint8_t>(std::clamp(new_pixel, 0.0, 255.0));
    }
  }

  return true;
}

bool BatkovFContrastEnhLinHistStretchSEQ::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace batkov_f_contrast_enh_lin_hist_stretch
