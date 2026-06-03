#include "romanova_v_linear_histogram_stretch/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "romanova_v_linear_histogram_stretch/common/include/common.hpp"

namespace romanova_v_linear_histogram_stretch_threads {

RomanovaVLinHistogramStretchSEQ::RomanovaVLinHistogramStretchSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool RomanovaVLinHistogramStretchSEQ::ValidationImpl() {
  return !GetInput().empty();
}

bool RomanovaVLinHistogramStretchSEQ::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return !GetOutput().empty();
}

bool RomanovaVLinHistogramStretchSEQ::RunImpl() {
  const InType &in = GetInput();
  OutType &out = GetOutput();

  auto [min_it, max_it] = std::ranges::minmax_element(in);
  uint8_t min_v = *min_it;
  uint8_t max_v = *max_it;

  if (min_v == max_v) {
    out = in;
    return true;
  }

  const uint8_t diff = max_v - min_v;

  for (size_t i = 0; i < in.size(); i++) {
    uint8_t pix = in[i];
    out[i] = (pix - min_v) / diff * 255;
  }

  return true;
}

bool RomanovaVLinHistogramStretchSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace romanova_v_linear_histogram_stretch_threads
