#include "batkov_f_contrast_enh_lin_hist_stretch/stl/include/ops_stl.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

#include "batkov_f_contrast_enh_lin_hist_stretch/common/include/common.hpp"
#include "util/include/util.hpp"

namespace batkov_f_contrast_enh_lin_hist_stretch {

namespace {

std::pair<uint8_t, uint8_t> FindMinMaxParallel(const InType &input, size_t num_threads) {
  const size_t n = input.size();
  const size_t block = n / num_threads;

  std::vector<uint8_t> mins(num_threads, std::numeric_limits<uint8_t>::max());
  std::vector<uint8_t> maxs(num_threads, std::numeric_limits<uint8_t>::min());

  {
    std::vector<std::thread> threads;
    threads.reserve(num_threads);
    for (size_t thread_index = 0; thread_index < num_threads; ++thread_index) {
      const size_t begin = thread_index * block;
      const size_t end = (thread_index == num_threads - 1) ? n : begin + block;

      threads.emplace_back([&, thread_index, begin, end]() {
        uint8_t local_min = std::numeric_limits<uint8_t>::max();
        uint8_t local_max = std::numeric_limits<uint8_t>::min();
        for (size_t i = begin; i < end; ++i) {
          local_min = std::min(local_min, input[i]);
          local_max = std::max(local_max, input[i]);
        }
        mins[thread_index] = local_min;
        maxs[thread_index] = local_max;
      });
    }

    for (auto &th : threads) {
      th.join();
    }
  }

  return {*std::ranges::min_element(mins), *std::ranges::max_element(maxs)};
}

std::pair<uint8_t, uint8_t> FindMinMax(const InType &input, size_t parallel_threshold, size_t num_threads) {
  if (input.size() < parallel_threshold || num_threads <= 1) {
    const auto [min_it, max_it] = std::ranges::minmax_element(input);
    return {*min_it, *max_it};
  }

  return FindMinMaxParallel(input, num_threads);
}

std::array<uint8_t, 256> BuildStretchLut(float a, float b) {
  std::array<uint8_t, 256> lut{};
  for (size_t pixel = 0; pixel < 256; ++pixel) {
    lut.at(pixel) = static_cast<uint8_t>(std::clamp((a * static_cast<float>(pixel)) + b, 0.0F, 255.0F));
  }
  return lut;
}

void ApplyStretchParallel(const InType &input, OutType &output, size_t num_threads,
                          const std::array<uint8_t, 256> &lut) {
  const size_t n = input.size();
  const size_t block = n / num_threads;

  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  for (size_t thread_index = 0; thread_index < num_threads; ++thread_index) {
    const size_t begin = thread_index * block;
    const size_t end = (thread_index == num_threads - 1) ? n : begin + block;

    threads.emplace_back([&, begin, end]() {
      for (size_t i = begin; i < end; ++i) {
        output[i] = lut.at(input[i]);
      }
    });
  }

  for (auto &th : threads) {
    th.join();
  }
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
  auto &input = GetInput();
  auto &output = GetOutput();

  constexpr size_t kParallelMinMaxThreshold = 100000;
  const size_t num_threads = static_cast<size_t>(std::max(1, ppc::util::GetNumThreads()));

  uint8_t min_el{};
  uint8_t max_el{};
  std::tie(min_el, max_el) = FindMinMax(input, kParallelMinMaxThreshold, num_threads);

  if (min_el == max_el) {
    std::ranges::copy(input, output.begin());
    return true;
  }

  const float a = 255.0F / static_cast<float>(max_el - min_el);
  const float b = -a * static_cast<float>(min_el);
  const auto lut = BuildStretchLut(a, b);
  ApplyStretchParallel(input, output, num_threads, lut);
  return true;
}

bool BatkovFContrastEnhLinHistStretchSTL::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace batkov_f_contrast_enh_lin_hist_stretch
