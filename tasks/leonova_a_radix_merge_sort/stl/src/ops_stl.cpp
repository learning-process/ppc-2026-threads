#include "leonova_a_radix_merge_sort/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <future>
#include <vector>

#include "leonova_a_radix_merge_sort/common/include/common.hpp"
#include "util/include/util.hpp"

namespace leonova_a_radix_merge_sort {

LeonovaARadixMergeSortSTL::LeonovaARadixMergeSortSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<int64_t>(GetInput().size());
}

bool LeonovaARadixMergeSortSTL::ValidationImpl() {
  return !GetInput().empty();
}

bool LeonovaARadixMergeSortSTL::PreProcessingImpl() {
  return true;
}

bool LeonovaARadixMergeSortSTL::RunImpl() {
  if (!ValidationImpl()) {
    return false;
  }
  GetOutput() = GetInput();
  if (GetOutput().size() > 1) {
    RadixMergeSort(GetOutput(), 0, GetOutput().size());
  }
  return true;
}

bool LeonovaARadixMergeSortSTL::PostProcessingImpl() {
  return true;
}

inline uint64_t LeonovaARadixMergeSortSTL::ToUnsignedValue(int64_t value) {
  return static_cast<uint64_t>(value) ^ kSignBitMask;
}

void LeonovaARadixMergeSortSTL::FillKeys(const std::vector<int64_t> &arr, size_t left, size_t size,
                                         std::vector<uint64_t> &keys) {
  for (size_t index = 0; index < size; ++index) {
    keys[index] = ToUnsignedValue(arr[left + index]);
  }
}

void LeonovaARadixMergeSortSTL::CountBytes(const std::vector<uint64_t> &keys, size_t begin, size_t end, int shift,
                                           std::vector<size_t> &counts) {
  std::ranges::fill(counts, 0);
  for (size_t index = begin; index < end; ++index) {
    ++counts[(keys[index] >> shift) & 0xFFU];
  }
}

void LeonovaARadixMergeSortSTL::ScatterBytes(const std::vector<uint64_t> &src_keys, const std::vector<int64_t> &src_arr,
                                             size_t left, size_t begin, size_t end, int shift,
                                             std::vector<size_t> &offsets, std::vector<int64_t> &dst_arr,
                                             std::vector<uint64_t> &dst_keys) {
  for (size_t index = begin; index < end; ++index) {
    const size_t bucket = (src_keys[index] >> shift) & 0xFFU;
    const size_t pos = offsets[bucket]++;
    dst_arr[pos] = src_arr[left + index];
    dst_keys[pos] = src_keys[index];
  }
}

void LeonovaARadixMergeSortSTL::RadixSort(std::vector<int64_t> &arr, size_t left, size_t right) {
  const size_t size = right - left;
  if (size <= 1) {
    return;
  }

  std::vector<uint64_t> keys(size);
  std::vector<int64_t> tmp_arr(size);
  std::vector<uint64_t> tmp_keys(size);
  std::vector<size_t> counts(kNumCounters);
  std::vector<size_t> offsets(kNumCounters);

  FillKeys(arr, left, size, keys);

  for (int byte_pos = 0; byte_pos < kNumBytes; ++byte_pos) {
    const int shift = byte_pos * kByteSize;

    CountBytes(keys, 0, size, shift, counts);
    size_t sum = 0;
    for (size_t bndex = 0; bndex < kNumCounters; ++bndex) {
      offsets[bndex] = sum;
      sum += counts[bndex];
    }

    ScatterBytes(keys, arr, left, 0, size, shift, offsets, tmp_arr, tmp_keys);

    std::ranges::copy(tmp_arr, arr.begin() + static_cast<std::ptrdiff_t>(left));
    keys.swap(tmp_keys);
  }
}

void LeonovaARadixMergeSortSTL::SimpleMerge(std::vector<int64_t> &arr, size_t left, size_t mid, size_t right) {
  std::vector<int64_t> merged(right - left);

  size_t i = left;
  size_t j = mid;
  size_t k = 0;

  while (i < mid && j < right) {
    merged[k++] = (arr[i] <= arr[j]) ? arr[i++] : arr[j++];
  }
  while (i < mid) {
    merged[k++] = arr[i++];
  }
  while (j < right) {
    merged[k++] = arr[j++];
  }

  std::ranges::copy(merged, arr.begin() + static_cast<std::ptrdiff_t>(left));
}

void LeonovaARadixMergeSortSTL::RadixMergeSort(std::vector<int64_t> &arr, size_t left, size_t right) {
  const size_t size = right - left;
  const auto hw = static_cast<size_t>(std::max(1, ppc::util::GetNumThreads()));
  const size_t num_threads = std::min(hw, size);
  std::vector<size_t> boundaries(num_threads + 1);
  const size_t chunk = (size + num_threads - 1) / num_threads;
  for (size_t tndex = 0; tndex <= num_threads; ++tndex) {
    boundaries[tndex] = left + std::min(tndex * chunk, size);
  }

  std::vector<std::future<void>> futures;
  futures.reserve(num_threads);

  for (size_t tndex = 0; tndex < num_threads; ++tndex) {
    const size_t l = boundaries[tndex];
    const size_t r = boundaries[tndex + 1];
    futures.push_back(std::async(std::launch::async, [&arr, l, r] { RadixSort(arr, l, r); }));
  }

  for (auto &f : futures) {
    f.get();
  }

  size_t step = chunk;
  while (step < size) {
    const size_t double_step = step * 2;
    size_t pos = left;
    while (pos < right) {
      const size_t mid = std::min(pos + step, right);
      const size_t end = std::min(pos + double_step, right);
      if (mid < end) {
        SimpleMerge(arr, pos, mid, end);
      }
      pos += double_step;
    }
    step = double_step;
  }
}

}  // namespace leonova_a_radix_merge_sort
