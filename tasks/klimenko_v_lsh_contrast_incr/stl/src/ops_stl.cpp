#include "klimenko_v_lsh_contrast_incr/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstddef>
#include <thread>
#include <vector>

#include "klimenko_v_lsh_contrast_incr/common/include/common.hpp"
#include "util/include/util.hpp"

namespace klimenko_v_lsh_contrast_incr {

namespace {

void GetThreadRange(size_t tid, size_t total, size_t num_t, size_t &b, size_t &e) {
  size_t chunk = total / num_t;
  b = tid * chunk;
  e = (tid == num_t - 1) ? total : b + chunk;
}

}  // namespace

KlimenkoVLSHContrastIncrSTL::KlimenkoVLSHContrastIncrSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool KlimenkoVLSHContrastIncrSTL::ValidationImpl() {
  return !GetInput().empty();
}

bool KlimenkoVLSHContrastIncrSTL::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return true;
}

bool KlimenkoVLSHContrastIncrSTL::RunImpl() {
  const auto &input = GetInput();
  auto &output = GetOutput();
  const size_t total_size = input.size();

  if (total_size == 0) {
    return false;
  }

  const int num_t = ppc::util::GetNumThreads();
  std::vector<std::thread> threads;
  threads.reserve(num_t);

  std::vector<int> local_min(num_t);
  std::vector<int> local_max(num_t);

  for (int tid = 0; tid < num_t; tid++) {
    threads.emplace_back([&, tid]() {
      size_t begin = 0;
      size_t end = 0;
      GetThreadRange(tid, total_size, num_t, begin, end);

      int current_min = input[begin];
      int current_max = input[begin];

      for (size_t i = begin; i < end; i++) {
        current_min = std::min(current_min, input[i]);
        current_max = std::max(current_max, input[i]);
      }

      local_min[tid] = current_min;
      local_max[tid] = current_max;
    });
  }

  for (auto &th : threads) {
    th.join();
  }
  threads.clear();

  const int global_min = *std::ranges::min_element(local_min);
  const int global_max = *std::ranges::max_element(local_max);

  if (global_max == global_min) {
    std::ranges::copy(input, output.begin());
    return true;
  }

  for (int tid = 0; tid < num_t; tid++) {
    threads.emplace_back([&, tid]() {
      size_t begin = 0;
      size_t end = 0;
      GetThreadRange(tid, total_size, num_t, begin, end);

      for (size_t i = begin; i < end; i++) {
        output[i] = ((input[i] - global_min) * 255) / (global_max - global_min);
      }
    });
  }

  for (auto &th : threads) {
    th.join();
  }

  return true;
}

bool KlimenkoVLSHContrastIncrSTL::PostProcessingImpl() {
  return GetOutput().size() == GetInput().size();
}

}  // namespace klimenko_v_lsh_contrast_incr
