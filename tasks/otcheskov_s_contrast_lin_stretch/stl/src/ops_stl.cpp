#include "otcheskov_s_contrast_lin_stretch/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <thread>
#include <vector>

#include "otcheskov_s_contrast_lin_stretch/common/include/common.hpp"
#include "util/include/util.hpp"

namespace otcheskov_s_contrast_lin_stretch {

OtcheskovSContrastLinStretchSTL::OtcheskovSContrastLinStretchSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool OtcheskovSContrastLinStretchSTL::ValidationImpl() {
  return !GetInput().empty();
}

bool OtcheskovSContrastLinStretchSTL::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return GetOutput().size() == GetInput().size();
}

bool OtcheskovSContrastLinStretchSTL::RunImpl() {
  if (GetInput().empty()) {
    return false;
  }

  const InType &input = GetInput();
  OutType &output = GetOutput();
  const size_t size = input.size();

  MinMax result = ComputeMinMax(input);
  if (result.min == result.max) {
    const size_t threshold_size = 1000000;
    if (size > threshold_size) {
      CopyInput(input, output);
    } else {
      for (size_t i = 0; i < size; ++i) {
        output[i] = input[i];
      }
    }
    return true;
  }
  const int min_i = static_cast<int>(result.min);
  const int range = static_cast<int>(result.max - min_i);
  LinearStretch(input, output, min_i, range);
  return true;
}

bool OtcheskovSContrastLinStretchSTL::PostProcessingImpl() {
  return true;
}

OtcheskovSContrastLinStretchSTL::MinMax OtcheskovSContrastLinStretchSTL::ComputeMinMax(const InType &input) {
  const size_t size = input.size();
  const size_t num_threads = std::min<size_t>(static_cast<size_t>(ppc::util::GetNumThreads()), size);

  const size_t block = size / num_threads;
  std::vector<MinMax> local(num_threads, {255, 0});

  {
    std::vector<std::jthread> threads;
    threads.reserve(num_threads);
    for (size_t tid = 0; tid < num_threads; ++tid) {
      size_t begin = tid * block;
      size_t end = (tid == num_threads - 1) ? size : begin + block;

      threads.emplace_back([&, tid, begin, end]() {
        auto [min, max] = std::ranges::minmax_element(input.begin() + static_cast<std::ptrdiff_t>(begin),
                                                      input.begin() + static_cast<std::ptrdiff_t>(end));
        local[tid] = {.min = *min, .max = *max};
      });
    }
  }

  MinMax result{.min = 255, .max = 0};
  for (size_t tid = 0; tid < num_threads; ++tid) {
    result.min = std::min(result.min, local[tid].min);
    result.max = std::max(result.max, local[tid].max);
  }
  return result;
}

void OtcheskovSContrastLinStretchSTL::CopyInput(const InType &input, OutType &output) {
  const size_t size = input.size();
  const size_t num_threads = std::min<size_t>(static_cast<size_t>(ppc::util::GetNumThreads()), size);
  const size_t block = size / num_threads;

  std::vector<std::jthread> threads;
  threads.reserve(num_threads);
  for (size_t tid = 0; tid < num_threads; ++tid) {
    size_t begin = tid * block;
    size_t end = (tid == num_threads - 1) ? size : begin + block;

    threads.emplace_back([&, begin, end]() {
      for (size_t i = begin; i < end; ++i) {
        output[i] = input[i];
      }
    });
  }
}

void OtcheskovSContrastLinStretchSTL::LinearStretch(const InType &input, OutType &output, int min_i, int range) {
  const size_t size = input.size();
  const size_t num_threads = std::min<size_t>(static_cast<size_t>(ppc::util::GetNumThreads()), size);
  const size_t block = size / num_threads;

  std::vector<std::jthread> threads;
  threads.reserve(num_threads);
  for (size_t tid = 0; tid < num_threads; ++tid) {
    size_t begin = tid * block;
    size_t end = (tid == num_threads - 1) ? size : begin + block;

    threads.emplace_back([&, begin, end, min_i, range]() {
      for (size_t i = begin; i < end; ++i) {
        int pixel = static_cast<int>(input[i]);
        int value = (pixel - min_i) * 255 / (range);
        value = std::clamp(value, 0, 255);
        output[i] = static_cast<uint8_t>(value);
      }
    });
  }
}

}  // namespace otcheskov_s_contrast_lin_stretch
