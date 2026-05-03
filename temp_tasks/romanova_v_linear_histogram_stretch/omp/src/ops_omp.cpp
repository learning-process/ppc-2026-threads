#include "romanova_v_linear_histogram_stretch/omp/include/ops_omp.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "romanova_v_linear_histogram_stretch/common/include/common.hpp"

namespace romanova_v_linear_histogram_stretch_threads {

RomanovaVLinHistogramStretchOMP::RomanovaVLinHistogramStretchOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool RomanovaVLinHistogramStretchOMP::ValidationImpl() {
  return !GetInput().empty();
}

bool RomanovaVLinHistogramStretchOMP::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return !GetOutput().empty();
}

bool RomanovaVLinHistogramStretchOMP::RunImpl() {
  const InType &in = GetInput();
  OutType &out = GetOutput();
  size_t size = in.size();

  uint8_t min_v = 255;
  uint8_t max_v = 0;

#pragma omp parallel for default(none) shared(in, size) reduction(min : min_v) reduction(max : max_v)
  for (size_t i = 0; i < size; i++) {
    min_v = std::min(min_v, in[i]);
    max_v = std::max(max_v, in[i]);
  }

  if (min_v == max_v) {
#pragma omp parallel for default(none) shared(out, in, size)
    for (size_t i = 0; i < size; i++) {
      out[i] = in[i];
    }
    return true;
  }

  const uint8_t diff = max_v - min_v;
  const double ratio = 255.0 / diff;

#pragma omp parallel for default(none) shared(in, out, min_v, ratio, size)
  for (size_t i = 0; i < size; i++) {
    uint8_t pix = in[i];
    out[i] = (std::clamp(static_cast<int>((pix - min_v) * ratio), 0, 255));
  }

  return true;
}

bool RomanovaVLinHistogramStretchOMP::PostProcessingImpl() {
  return true;
}

}  // namespace romanova_v_linear_histogram_stretch_threads
