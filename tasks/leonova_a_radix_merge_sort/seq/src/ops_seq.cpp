#include "leonova_a_radix_merge_sort/seq/include/ops_seq.hpp"

#include <algorithm>  // для std::ranges::copy
#include <cstddef>    // для size_t
#include <cstdint>    // для uint64_t, int64_t, uint8_t
#include <vector>

#include "leonova_a_radix_merge_sort/common/include/common.hpp"

namespace leonova_a_radix_merge_sort {

LeonovaARadixMergeSortSEQ::LeonovaARadixMergeSortSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<int64_t>(GetInput().size());
}

bool LeonovaARadixMergeSortSEQ::ValidationImpl() {
  return !GetInput().empty();
}

bool LeonovaARadixMergeSortSEQ::PreProcessingImpl() {
  return true;
}

bool LeonovaARadixMergeSortSEQ::RunImpl() {
  if (!ValidationImpl()) {
    return false;
  }
  GetOutput() = GetInput();
  if (GetOutput().size() > 1) {
    RadixMergeSort(GetOutput(), 0, GetOutput().size());
  }
  return true;
}

bool LeonovaARadixMergeSortSEQ::PostProcessingImpl() {
  return true;
}

uint64_t LeonovaARadixMergeSortSEQ::ToUnsigned(int64_t value) {
  auto unsigned_value = static_cast<uint64_t>(value);
  constexpr uint64_t kSignBitMask = 0x8000000000000000ULL;
  return unsigned_value ^ kSignBitMask;
}

void LeonovaARadixMergeSortSEQ::RadixSort(std::vector<int64_t> &arr, size_t left, size_t right) {
  size_t size = right - left;
  if (size <= 1) {
    return;
  }

  std::vector<uint64_t> keys(size);
  for (size_t index = 0; index < size; ++index) {
    keys[index] = ToUnsigned(arr[left + index]);
  }

  std::vector<int64_t> temp_arr(size);
  std::vector<uint64_t> temp_keys(size);

  for (int byte_pos = 0; byte_pos < kNumBytes; ++byte_pos) {
    std::vector<size_t> count(kNumCounters, 0);

    for (size_t index = 0; index < size; ++index) {
      uint8_t byte_val = (keys[index] >> (byte_pos * kByteSize)) & 0xFF;
      ++count[byte_val];
    }

    size_t total = 0;
    for (size_t &elem : count) {
      size_t old_count = elem;
      elem = total;
      total += old_count;
    }

    for (size_t index = 0; index < size; ++index) {
      uint8_t byte_val = (keys[index] >> (byte_pos * kByteSize)) & 0xFF;
      size_t pos = count[byte_val]++;
      temp_arr[pos] = arr[left + index];
      temp_keys[pos] = keys[index];
    }

    auto left_it = arr.begin() + static_cast<typename std::vector<int64_t>::difference_type>(left);
    std::ranges::copy(temp_arr, left_it);
    keys.swap(temp_keys);
  }
}

void LeonovaARadixMergeSortSEQ::SimpleMerge(std::vector<int64_t> &arr, 
                                            size_t left, size_t mid, size_t right) {
  size_t left_size = mid - left;
  size_t right_size = right - mid;

  std::vector<int64_t> merged(left_size + right_size);

  size_t i = 0;
  size_t j = 0;
  size_t k = 0;

  while (i < left_size && j < right_size) {
    if (arr[left + i] <= arr[mid + j]) {
      merged[k++] = arr[left + i++];
    } else {
      merged[k++] = arr[mid + j++];
    }
  }

  while (i < left_size) {
    merged[k++] = arr[left + i++];
  }
  while (j < right_size) {
    merged[k++] = arr[mid + j++];
  }

  auto left_it = arr.begin() + static_cast<typename std::vector<int64_t>::difference_type>(left);
  std::ranges::copy(merged, left_it);
}

void LeonovaARadixMergeSortSEQ::RadixMergeSort(std::vector<int64_t> &arr, 
                                               size_t left, size_t right) {
  struct SortTask {
    size_t left;
    size_t right;
    bool sorted;
  };

  std::vector<SortTask> stack;
  stack.reserve(128);
  stack.push_back({left, right, false});

  while (!stack.empty()) {
    SortTask current = stack.back();
    stack.pop_back();

    size_t size = current.right - current.left;

    if (size <= 1) {
      continue;
    }

    if (size <= kRadixThreshold) {
      RadixSort(arr, current.left, current.right);
      continue;
    }

    if (!current.sorted) {
      size_t mid = current.left + (size / 2);

      stack.push_back({current.left, current.right, true});
      stack.push_back({mid, current.right, false});
      stack.push_back({current.left, mid, false});
    } else {
      size_t mid = current.left + (size / 2);
      SimpleMerge(arr, current.left, mid, current.right);
    }
  }
}
}  // namespace leonova_a_radix_merge_sort