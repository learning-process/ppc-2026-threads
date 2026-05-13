#include "romanova_v_linear_histogram_stretch/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <thread>
#include <vector>

#include "romanova_v_linear_histogram_stretch/common/include/common.hpp"
#include "util/include/util.hpp"

namespace romanova_v_linear_histogram_stretch_threads {

void RomanovaVLinHistogramStretchSTL::GetThreadRange(size_t thid, size_t total, size_t num_th, size_t &beg,
                                                     size_t &en) {
  size_t chunk = total / num_th;
  beg = thid * chunk;
  en = (thid == num_th - 1) ? total : (beg + chunk);
}

void RomanovaVLinHistogramStretchSTL::FindMinMax(size_t begin, size_t end, uint8_t &curr_min, uint8_t &curr_max,
                                                 const InType &in) {
  for (size_t i = begin; i < end; ++i) {
    uint8_t val = in[i];
    curr_min = std::min(curr_min, val);
    curr_max = std::max(curr_max, val);
  }
}

RomanovaVLinHistogramStretchSTL::RomanovaVLinHistogramStretchSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool RomanovaVLinHistogramStretchSTL::ValidationImpl() {
  return !GetInput().empty();
}

bool RomanovaVLinHistogramStretchSTL::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return !GetOutput().empty();
}

bool RomanovaVLinHistogramStretchSTL::RunImpl() {
  const InType &in = GetInput();
  OutType &out = GetOutput();
  size_t size = in.size();

  const int num_th = ppc::util::GetNumThreads();
  std::vector<std::thread> threads;
  threads.reserve(num_th);

  std::vector<uint8_t> local_min(num_th, 255);
  std::vector<uint8_t> local_max(num_th, 0);

  for (int thid = 0; thid < num_th; thid++) {
    threads.emplace_back([&, thid]() {
      size_t begin = 0;
      size_t end = 0;

      GetThreadRange(thid, size, num_th, begin, end);

      uint8_t current_min = 255;
      uint8_t current_max = 0;

      FindMinMax(begin, end, current_min, current_max, in);

      local_min[thid] = current_min;
      local_max[thid] = current_max;
    });
  }

  for (auto &th : threads) {
    th.join();
  }
  threads.clear();

  uint8_t min_v = local_min[0];
  uint8_t max_v = local_max[0];
  for (int thid = 1; thid < num_th; thid++) {
    min_v = std::min(min_v, local_min[thid]);
    max_v = std::max(max_v, local_max[thid]);
  }

  if (min_v == max_v) {
    std::ranges::copy(in, out.begin());
    return true;
  }

  const uint8_t diff = max_v - min_v;
  const double ratio = 255.0 / diff;

  for (int thid = 0; thid < num_th; thid++) {
    threads.emplace_back([&, thid]() {
      size_t begin = 0;
      size_t end = 0;

      GetThreadRange(thid, size, num_th, begin, end);

      for (size_t i = begin; i < end; i++) {
        out[i] = (std::clamp(static_cast<int>((in[i] - min_v) * ratio), 0, 255));
      }
    });
  }

  for (auto &th : threads) {
    th.join();
  }
  threads.clear();

  return true;
}

bool RomanovaVLinHistogramStretchSTL::PostProcessingImpl() {
  return true;
}

}  // namespace romanova_v_linear_histogram_stretch_threads
