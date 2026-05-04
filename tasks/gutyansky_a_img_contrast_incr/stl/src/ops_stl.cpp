#include "gutyansky_a_img_contrast_incr/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <thread>
#include <utility>
#include <vector>

#include "gutyansky_a_img_contrast_incr/common/include/common.hpp"

namespace gutyansky_a_img_contrast_incr {

namespace {
void WaitAll(std::vector<std::thread> &threads) {
  for (auto &thread : threads) {
    thread.join();
  }
}

std::pair<uint8_t, uint8_t> ComputeBounds(const std::vector<uint8_t> &input) {
  const size_t sz = input.size();
  const uint32_t num_threads = std::max(1U, std::thread::hardware_concurrency());
  const size_t chunk_sz = (sz + num_threads - 1) / num_threads;

  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  std::vector<std::pair<uint8_t, uint8_t>> local_res(
      num_threads, std::make_pair(static_cast<uint8_t>(255), static_cast<uint8_t>(0)));

  for (uint32_t tid = 0; tid < num_threads; tid++) {
    size_t from = tid * chunk_sz;
    size_t to = std::min(sz, from + chunk_sz);

    threads.emplace_back([&local_res, &input, from, to, tid]() {
      uint8_t local_min = 255;
      uint8_t local_max = 0;

      for (size_t i = from; i < to; ++i) {
        local_min = std::min(local_min, input[i]);
        local_max = std::max(local_max, input[i]);
      }

      local_res[tid] = std::make_pair(local_min, local_max);
    });
  }

  WaitAll(threads);

  uint8_t lower_bound = 255;
  uint8_t upper_bound = 0;
  for (const auto &result : local_res) {
    lower_bound = std::min(lower_bound, result.first);
    upper_bound = std::max(upper_bound, result.second);
  }

  return std::make_pair(lower_bound, upper_bound);
}

void CopyMultiThread(const std::vector<uint8_t> &input, std::vector<uint8_t> &output) {
  const size_t sz = input.size();
  const uint32_t num_threads = std::max(1U, std::thread::hardware_concurrency());
  const size_t chunk_sz = (sz + num_threads - 1) / num_threads;

  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  for (uint32_t tid = 0; tid < num_threads; tid++) {
    size_t from = tid * chunk_sz;
    size_t to = std::min(from + chunk_sz, sz);

    threads.emplace_back([&output, &input, from, to]() {
      for (size_t idx = from; idx < to; ++idx) {
        output[idx] = input[idx];
      }
    });
  }

  WaitAll(threads);
}

void TransformMultiThread(const std::vector<uint8_t> &input, std::vector<uint8_t> &output, uint8_t lower_bound,
                          uint8_t delta) {
  const size_t sz = input.size();
  const uint32_t num_threads = std::max(1U, std::thread::hardware_concurrency());
  const size_t chunk_sz = (sz + num_threads - 1) / num_threads;

  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  constexpr uint16_t kMaxUint8 = std::numeric_limits<uint8_t>::max();

  for (unsigned int tid = 0; tid < num_threads; tid++) {
    size_t from = tid * chunk_sz;
    size_t to = std::min(from + chunk_sz, sz);

    threads.emplace_back([&output, &input, from, to, lower_bound, delta]() {
      for (size_t idx = from; idx < to; ++idx) {
        uint16_t old_value = input[idx];
        uint16_t new_value = (kMaxUint8 * (old_value - lower_bound)) / delta;
        output[idx] = static_cast<uint8_t>(new_value);
      }
    });
  }

  WaitAll(threads);
}

}  // namespace

GutyanskyAImgContrastIncrSTL::GutyanskyAImgContrastIncrSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool GutyanskyAImgContrastIncrSTL::ValidationImpl() {
  return !GetInput().empty();
}

bool GutyanskyAImgContrastIncrSTL::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return true;
}

bool GutyanskyAImgContrastIncrSTL::RunImpl() {
  const auto &input = GetInput();
  auto &output = GetOutput();

  const auto bounds = ComputeBounds(input);
  const uint8_t delta = bounds.second - bounds.first;

  if (delta == 0) {
    CopyMultiThread(input, output);
  } else {
    TransformMultiThread(input, output, bounds.first, delta);
  }

  return true;
}

bool GutyanskyAImgContrastIncrSTL::PostProcessingImpl() {
  return true;
}

}  // namespace gutyansky_a_img_contrast_incr
