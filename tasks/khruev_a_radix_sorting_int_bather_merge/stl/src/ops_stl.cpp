#include "khruev_a_radix_sorting_int_bather_merge/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <thread>
#include <vector>

#include "khruev_a_radix_sorting_int_bather_merge/common/include/common.hpp"

namespace khruev_a_radix_sorting_int_bather_merge {

// Анонимное пространство имён для вспомогательных функций, чтобы они не засоряли глобальный скоуп
namespace {

template <typename Func>
void RunInThreads(size_t total_tasks, size_t max_threads, Func task_func) {
  size_t num_threads = std::min(max_threads, total_tasks);
  if (num_threads <= 1) {
    for (size_t task_idx = 0; task_idx < total_tasks; ++task_idx) {
      task_func(task_idx);
    }
    return;
  }

  std::vector<std::thread> threads;
  threads.reserve(num_threads);
  size_t chunk_size = (total_tasks + num_threads - 1) / num_threads;

  for (size_t thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    size_t start_idx = thread_idx * chunk_size;
    size_t end_idx = std::min(start_idx + chunk_size, total_tasks);

    if (start_idx < end_idx) {
      threads.emplace_back([start_idx, end_idx, &task_func]() {
        for (size_t task_idx = start_idx; task_idx < end_idx; ++task_idx) {
          task_func(task_idx);
        }
      });
    }
  }

  for (auto &thread_obj : threads) {
    thread_obj.join();
  }
}

}  // namespace

void KhruevARadixSortingIntBatherMergeSTL::CompareExchange(std::vector<int> &a, size_t i, size_t j) {
  if (a[i] > a[j]) {
    std::swap(a[i], a[j]);
  }
}

void KhruevARadixSortingIntBatherMergeSTL::RadixSort(std::vector<int> &arr) {
  const int bits = 8;
  const int buckets = 1 << bits;
  const int mask = buckets - 1;
  const int passes = 32 / bits;

  std::vector<int> temp(arr.size());
  std::vector<int> *src = &arr;
  std::vector<int> *dst = &temp;

  for (int pass = 0; pass < passes; pass++) {
    std::vector<int> count(buckets, 0);
    int shift = pass * bits;

    for (int x : *src) {
      uint32_t ux = static_cast<uint32_t>(x) ^ 0x80000000U;
      uint32_t digit = (ux >> shift) & mask;
      count[digit]++;
    }

    for (int i = 1; i < buckets; i++) {
      count[i] += count[i - 1];
    }

    for (size_t i = src->size(); i-- > 0;) {
      uint32_t ux = static_cast<uint32_t>((*src)[i]) ^ 0x80000000U;
      uint32_t digit = (ux >> shift) & mask;
      (*dst)[--count[digit]] = (*src)[i];
    }

    std::swap(src, dst);
  }
}

void KhruevARadixSortingIntBatherMergeSTL::OddEvenMerge(std::vector<int> &a, size_t n) {
  if (n < 2) {
    return;
  }

  size_t max_threads = std::thread::hardware_concurrency();
  if (max_threads == 0) {
    max_threads = 4;
  }

  size_t po = n / 2;

  RunInThreads(po, max_threads, [&](size_t i) { CompareExchange(a, i, i + po); });

  po >>= 1;

  for (; po > 0; po >>= 1) {
    size_t num_blocks = (n - (2 * po)) / (2 * po);

    RunInThreads(num_blocks, max_threads, [&](size_t block_idx) {
      size_t i = po + (block_idx * 2 * po);
      for (size_t j = 0; j < po; ++j) {
        CompareExchange(a, i + j, i + j + po);
      }
    });
  }
}

KhruevARadixSortingIntBatherMergeSTL::KhruevARadixSortingIntBatherMergeSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool KhruevARadixSortingIntBatherMergeSTL::ValidationImpl() {
  return !GetInput().empty();
}

bool KhruevARadixSortingIntBatherMergeSTL::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return true;
}

bool KhruevARadixSortingIntBatherMergeSTL::RunImpl() {
  std::vector<int> data = GetInput();
  size_t original_size = data.size();

  if (original_size <= 1) {
    GetOutput() = data;
    return true;
  }

  size_t pow2 = 1;
  while (pow2 < original_size) {
    pow2 <<= 1;
  }

  data.resize(pow2, std::numeric_limits<int>::max());

  size_t half = pow2 / 2;
  auto half_dist = static_cast<std::ptrdiff_t>(half);

  std::vector<int> left(data.begin(), data.begin() + half_dist);
  std::vector<int> right(data.begin() + half_dist, data.end());

  std::thread left_thread([&]() { RadixSort(left); });

  RadixSort(right);

  left_thread.join();

  std::ranges::copy(left, data.begin());
  std::ranges::copy(right, data.begin() + half_dist);

  OddEvenMerge(data, data.size());

  data.resize(original_size);
  GetOutput() = data;

  return true;
}

bool KhruevARadixSortingIntBatherMergeSTL::PostProcessingImpl() {
  return GetOutput().size() == GetInput().size();
}

}  // namespace khruev_a_radix_sorting_int_bather_merge
