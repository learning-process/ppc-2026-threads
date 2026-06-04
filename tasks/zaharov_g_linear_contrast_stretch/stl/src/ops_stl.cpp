#include "zaharov_g_linear_contrast_stretch/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <thread>
#include <vector>

#include "util/include/util.hpp"
#include "zaharov_g_linear_contrast_stretch/common/include/common.hpp"

namespace zaharov_g_linear_contrast_stretch {

namespace {

struct MinMax {
  int min;
  int max;
};

std::size_t GetThreadCount(std::size_t size) {
  const auto requested_threads = static_cast<std::size_t>(std::max(1, ppc::util::GetNumThreads()));
  return std::min(size, requested_threads);
}

MinMax FindMinMax(const InType &input) {
  const std::size_t thread_count = GetThreadCount(input.size());
  std::vector<MinMax> local_minmax(thread_count,
                                   {std::numeric_limits<uint8_t>::max(), std::numeric_limits<uint8_t>::min()});
  std::vector<std::thread> threads;
  threads.reserve(thread_count);

  for (std::size_t thread_index = 0; thread_index < thread_count; ++thread_index) {
    const std::size_t begin = input.size() * thread_index / thread_count;
    const std::size_t end = input.size() * (thread_index + 1) / thread_count;
    threads.emplace_back([&input, &local_minmax, begin, end, thread_index]() {
      MinMax current{.min = std::numeric_limits<uint8_t>::max(), .max = std::numeric_limits<uint8_t>::min()};
      for (std::size_t i = begin; i < end; ++i) {
        const int value = static_cast<int>(input[i]);
        current.min = std::min(current.min, value);
        current.max = std::max(current.max, value);
      }
      local_minmax[thread_index] = current;
    });
  }

  for (auto &thread : threads) {
    thread.join();
  }

  MinMax result{.min = std::numeric_limits<uint8_t>::max(), .max = std::numeric_limits<uint8_t>::min()};
  for (const auto &current : local_minmax) {
    result.min = std::min(result.min, current.min);
    result.max = std::max(result.max, current.max);
  }
  return result;
}

void StretchRange(const InType &input, OutType &output, std::size_t begin, std::size_t end, int min_el, int max_el) {
  if (max_el > min_el) {
    const int denom = max_el - min_el;
    for (std::size_t i = begin; i < end; ++i) {
      const int value = (static_cast<int>(input[i]) - min_el) * std::numeric_limits<uint8_t>::max() / denom;
      output[i] = static_cast<uint8_t>(std::clamp(value, 0, 255));
    }
  } else {
    std::copy(input.begin() + static_cast<std::ptrdiff_t>(begin), input.begin() + static_cast<std::ptrdiff_t>(end),
              output.begin() + static_cast<std::ptrdiff_t>(begin));
  }
}

}  // namespace

ZaharovGLinContrStrSTL::ZaharovGLinContrStrSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool ZaharovGLinContrStrSTL::ValidationImpl() {
  return !GetInput().empty();
}

bool ZaharovGLinContrStrSTL::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return true;
}

bool ZaharovGLinContrStrSTL::RunImpl() {
  const InType &input = GetInput();
  OutType &output = GetOutput();

  const MinMax minmax = FindMinMax(input);
  const std::size_t thread_count = GetThreadCount(input.size());
  std::vector<std::thread> threads;
  threads.reserve(thread_count);

  for (std::size_t thread_index = 0; thread_index < thread_count; ++thread_index) {
    const std::size_t begin = input.size() * thread_index / thread_count;
    const std::size_t end = input.size() * (thread_index + 1) / thread_count;
    threads.emplace_back(
        [&input, &output, begin, end, minmax]() { StretchRange(input, output, begin, end, minmax.min, minmax.max); });
  }

  for (auto &thread : threads) {
    thread.join();
  }

  return true;
}

bool ZaharovGLinContrStrSTL::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace zaharov_g_linear_contrast_stretch
