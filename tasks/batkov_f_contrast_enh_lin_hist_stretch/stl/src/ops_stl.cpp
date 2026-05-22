#include "batkov_f_contrast_enh_lin_hist_stretch/stl/include/ops_stl.hpp"

#include <algorithm>
#include <array>
#include <barrier>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <thread>
#include <vector>

#include "batkov_f_contrast_enh_lin_hist_stretch/common/include/common.hpp"
#include "util/include/util.hpp"

namespace batkov_f_contrast_enh_lin_hist_stretch {

namespace {

constexpr size_t kParallelThreshold = 100000;
constexpr size_t kMinGrainSize = 1U << 16U;

struct Range {
  size_t begin;
  size_t end;

  Range(size_t begin_value, size_t end_value) noexcept : begin(begin_value), end(end_value) {}
};

struct MinMax {
  uint8_t min_val;
  uint8_t max_val;

  MinMax() noexcept : min_val(std::numeric_limits<uint8_t>::max()), max_val(std::numeric_limits<uint8_t>::min()) {}

  MinMax(uint8_t min_value, uint8_t max_value) noexcept : min_val(min_value), max_val(max_value) {}
};

size_t ChooseThreadCount(size_t size, size_t requested_threads) noexcept {
  if (size < kParallelThreshold || requested_threads <= 1) {
    return 1;
  }

  const size_t by_grain = std::max<size_t>(1, size / kMinGrainSize);
  return std::max<size_t>(1, std::min({requested_threads, size, by_grain}));
}

Range MakeRange(size_t size, size_t thread_index, size_t thread_count) noexcept {
  return {(size * thread_index) / thread_count, (size * (thread_index + 1)) / thread_count};
}

std::array<uint8_t, 256> BuildStretchLut(uint8_t min_el, uint8_t max_el) {
  std::array<uint8_t, 256> lut{};

  const int min_value = static_cast<int>(min_el);
  const int max_value = static_cast<int>(max_el);
  const int range = max_value - min_value;

  for (size_t pixel = 0; pixel < lut.size(); ++pixel) {
    const int value = static_cast<int>(pixel);

    if (value <= min_value) {
      lut.at(pixel) = 0;
    } else if (value >= max_value) {
      lut.at(pixel) = std::numeric_limits<uint8_t>::max();
    } else {
      lut.at(pixel) = static_cast<uint8_t>(((value - min_value) * 255) / range);
    }
  }

  return lut;
}

uint8_t GetLutValue(const std::array<uint8_t, 256> &lut, uint8_t pixel) noexcept {
  return *std::next(lut.cbegin(), static_cast<std::ptrdiff_t>(pixel));
}

MinMax FindLocalMinMax(const InType &input, const Range &range) {
  auto input_it = std::next(input.cbegin(), static_cast<std::ptrdiff_t>(range.begin));
  const auto input_last = std::next(input.cbegin(), static_cast<std::ptrdiff_t>(range.end));

  MinMax result;
  for (; input_it != input_last; ++input_it) {
    result.min_val = std::min(result.min_val, *input_it);
    result.max_val = std::max(result.max_val, *input_it);
  }

  return result;
}

MinMax CombineLocalMinMax(const std::vector<MinMax> &locals) noexcept {
  MinMax result;

  for (const auto &local : locals) {
    result.min_val = std::min(result.min_val, local.min_val);
    result.max_val = std::max(result.max_val, local.max_val);
  }

  return result;
}

void CopyRange(const InType &input, OutType &output, const Range &range) {
  const auto input_first = std::next(input.cbegin(), static_cast<std::ptrdiff_t>(range.begin));
  const auto input_last = std::next(input.cbegin(), static_cast<std::ptrdiff_t>(range.end));
  const auto output_first = std::next(output.begin(), static_cast<std::ptrdiff_t>(range.begin));

  std::copy(input_first, input_last, output_first);
}

void ApplyLutRange(const InType &input, OutType &output, const Range &range, const std::array<uint8_t, 256> &lut) {
  auto input_it = std::next(input.cbegin(), static_cast<std::ptrdiff_t>(range.begin));
  const auto input_last = std::next(input.cbegin(), static_cast<std::ptrdiff_t>(range.end));
  auto output_it = std::next(output.begin(), static_cast<std::ptrdiff_t>(range.begin));

  for (; input_it != input_last; ++input_it, ++output_it) {
    *output_it = GetLutValue(lut, *input_it);
  }
}

void ApplyResultRange(const InType &input, OutType &output, const Range &range, bool copy_only,
                      const std::array<uint8_t, 256> &lut) {
  if (copy_only) {
    CopyRange(input, output, range);
    return;
  }

  ApplyLutRange(input, output, range, lut);
}

void StretchSequential(const InType &input, OutType &output) {
  const auto minmax = std::ranges::minmax_element(input.cbegin(), input.cend());
  const uint8_t min_el = *minmax.max;
  const uint8_t max_el = *minmax.min;
  const Range full_range(0, input.size());

  if (min_el == max_el) {
    CopyRange(input, output, full_range);
    return;
  }

  const auto lut = BuildStretchLut(min_el, max_el);
  ApplyLutRange(input, output, full_range, lut);
}

void StretchParallel(const InType &input, OutType &output, size_t thread_count) {
  const size_t size = input.size();
  std::vector<MinMax> locals(thread_count);
  std::array<uint8_t, 256> lut{};
  bool copy_only = false;

  auto completion = [&locals, &lut, &copy_only]() noexcept {
    const MinMax global = CombineLocalMinMax(locals);
    copy_only = global.min_val == global.max_val;

    if (!copy_only) {
      lut = BuildStretchLut(global.min_val, global.max_val);
    }
  };

  std::barrier<decltype(completion)> barrier(static_cast<std::ptrdiff_t>(thread_count), completion);

  std::vector<std::thread> threads;
  threads.reserve(thread_count);

  for (size_t thread_index = 0; thread_index < thread_count; ++thread_index) {
    threads.emplace_back([&, thread_index]() {
      const Range range = MakeRange(size, thread_index, thread_count);
      locals.at(thread_index) = FindLocalMinMax(input, range);
      barrier.arrive_and_wait();
      ApplyResultRange(input, output, range, copy_only, lut);
    });
  }

  for (auto &thread : threads) {
    thread.join();
  }
}

void StretchContrast(const InType &input, OutType &output, size_t requested_threads) {
  const size_t thread_count = ChooseThreadCount(input.size(), requested_threads);

  if (thread_count == 1) {
    StretchSequential(input, output);
    return;
  }

  StretchParallel(input, output, thread_count);
}

}  // namespace

BatkovFContrastEnhLinHistStretchSTL::BatkovFContrastEnhLinHistStretchSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool BatkovFContrastEnhLinHistStretchSTL::ValidationImpl() {
  return !GetInput().empty();
}

bool BatkovFContrastEnhLinHistStretchSTL::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return true;
}

bool BatkovFContrastEnhLinHistStretchSTL::RunImpl() {
  const auto &input = GetInput();
  auto &output = GetOutput();

  const size_t requested_threads = static_cast<size_t>(std::max(1, ppc::util::GetNumThreads()));
  StretchContrast(input, output, requested_threads);

  return true;
}

bool BatkovFContrastEnhLinHistStretchSTL::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace batkov_f_contrast_enh_lin_hist_stretch
