#include "romanova_v_linear_histogram_stretch/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/tbb.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "romanova_v_linear_histogram_stretch/common/include/common.hpp"

namespace romanova_v_linear_histogram_stretch_threads {

RomanovaVLinHistogramStretchTBB::RomanovaVLinHistogramStretchTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool RomanovaVLinHistogramStretchTBB::ValidationImpl() {
  return !GetInput().empty();
}

bool RomanovaVLinHistogramStretchTBB::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return !GetOutput().empty();
}

bool RomanovaVLinHistogramStretchTBB::RunImpl() {
  const InType &in = GetInput();
  OutType &out = GetOutput();

  struct MM {
    uint8_t min;
    uint8_t max;
  };

  MM minmax = tbb::parallel_reduce(tbb::blocked_range<size_t>(0, in.size()), MM{.min = 255, .max = 0},
                                   [&](const tbb::blocked_range<size_t> &range, MM init) {
    for (size_t i = range.begin(); i != range.end(); i++) {
      init = MM{.min = std::min(init.min, in[i]), .max = std::max(init.max, in[i])};
    }
    return init;
  }, [](MM first, MM second) {
    return MM{.min = std::min(first.min, second.min), .max = std::max(first.max, second.max)};
  });

  if (minmax.min == minmax.max) {
    tbb::parallel_for(tbb::blocked_range<size_t>(0, in.size()), [&](const tbb::blocked_range<size_t> &range) {
      for (size_t i = range.begin(); i != range.end(); i++) {
        out[i] = in[i];
      }
    });
    return true;
  }

  const uint8_t diff = minmax.max - minmax.min;
  const double ratio = 255.0 / diff;

  tbb::parallel_for(tbb::blocked_range<size_t>(0, in.size()), [&](const tbb::blocked_range<size_t> &range) {
    for (size_t i = range.begin(); i != range.end(); i++) {
      out[i] = (std::clamp(static_cast<int>((in[i] - minmax.min) * ratio), 0, 255));
    }
  });

  return true;
}

bool RomanovaVLinHistogramStretchTBB::PostProcessingImpl() {
  return true;
}

}  // namespace romanova_v_linear_histogram_stretch_threads
