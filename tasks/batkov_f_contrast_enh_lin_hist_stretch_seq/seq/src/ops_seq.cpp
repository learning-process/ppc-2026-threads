#include "batkov_f_contrast_enh_lin_hist_stretch_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "batkov_f_contrast_enh_lin_hist_stretch_seq/common/include/common.hpp"

namespace batkov_f_contrast_enh_lin_hist_stretch_seq {

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
  auto [min_it, max_it] = std::ranges::minmax_element(GetInput());
  uint8_t min_el = *min_it;
  uint8_t max_el = *max_it;

  if (max_el > min_el) {
    double scale = 255.0 / (max_el - min_el);
    for (size_t i = 0; i < GetInput().size(); ++i) {
      GetOutput()[i] = static_cast<uint8_t>(std::min(255.0, (GetInput()[i] - min_el) * scale));
    }
  }

  return true;
}

bool BatkovFContrastEnhLinHistStretchSEQ::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace batkov_f_contrast_enh_lin_hist_stretch_seq
