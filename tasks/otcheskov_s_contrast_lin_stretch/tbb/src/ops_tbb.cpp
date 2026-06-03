#include "otcheskov_s_contrast_lin_stretch/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <util/include/util.hpp>
#include <vector>

#include "otcheskov_s_contrast_lin_stretch/common/include/common.hpp"

namespace otcheskov_s_contrast_lin_stretch {

OtcheskovSContrastLinStretchTBB::OtcheskovSContrastLinStretchTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool OtcheskovSContrastLinStretchTBB::ValidationImpl() {
  return !GetInput().empty();
}

bool OtcheskovSContrastLinStretchTBB::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return GetOutput().size() == GetInput().size();
}

bool OtcheskovSContrastLinStretchTBB::RunImpl() {
  if (GetInput().empty()) {
    return false;
  }

  const InType &input = GetInput();
  OutType &output = GetOutput();
  const size_t size = input.size();

  tbb::task_arena arena(ppc::util::GetNumThreads());

  MinMax result = ComputeMinMax(input, arena);
  if (result.min == result.max) {
    const size_t threshold_size = 1000000;
    if (size > threshold_size) {
      CopyInput(input, output, arena);
    } else {
      for (size_t i = 0; i < size; ++i) {
        output[i] = input[i];
      }
    }
    return true;
  }
  const int min_i = static_cast<int>(result.min);
  const int range = static_cast<int>(result.max - min_i);
  LinearStretch(input, output, min_i, range, arena);
  return true;
}

bool OtcheskovSContrastLinStretchTBB::PostProcessingImpl() {
  return true;
}

OtcheskovSContrastLinStretchTBB::MinMax OtcheskovSContrastLinStretchTBB::ComputeMinMax(const InType &input,
                                                                                       tbb::task_arena &arena) {
  MinMax result{};
  arena.execute([&] {
    result = tbb::parallel_reduce(tbb::blocked_range<size_t>(0, input.size()), MinMax{},
                                  [&](const tbb::blocked_range<size_t> &r, MinMax local) {
      for (size_t i = r.begin(); i != r.end(); ++i) {
        local.min = std::min(local.min, input[i]);
        local.max = std::max(local.max, input[i]);
      }
      return local;
    }, [](const MinMax &a, const MinMax &b) -> MinMax {
      return MinMax{.min = std::min(a.min, b.min), .max = std::max(a.max, b.max)};
    }, tbb::static_partitioner{});
  });
  return result;
}

void OtcheskovSContrastLinStretchTBB::CopyInput(const InType &input, OutType &output, tbb::task_arena &arena) {
  arena.execute([&] {
    tbb::parallel_for(tbb::blocked_range<size_t>(0, input.size()), [&](const tbb::blocked_range<size_t> &r) {
      for (size_t i = r.begin(); i != r.end(); ++i) {
        output[i] = input[i];
      }
    }, tbb::static_partitioner{});
  });
}

void OtcheskovSContrastLinStretchTBB::LinearStretch(const InType &input, OutType &output, int min_i, int range,
                                                    tbb::task_arena &arena) {
  arena.execute([&] {
    tbb::parallel_for(tbb::blocked_range<size_t>(0, input.size()), [&](const tbb::blocked_range<size_t> &r) {
      for (size_t i = r.begin(); i != r.end(); ++i) {
        int value = (static_cast<int>(input[i]) - min_i) * 255 / range;
        output[i] = static_cast<uint8_t>(std::clamp(value, 0, 255));
      }
    }, tbb::static_partitioner{});
  });
}

}  // namespace otcheskov_s_contrast_lin_stretch
