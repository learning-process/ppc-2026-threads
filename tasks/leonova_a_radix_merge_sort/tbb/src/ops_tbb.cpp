#include "leonova_a_radix_merge_sort/tbb/include/ops_tbb.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "leonova_a_radix_merge_sort/common/include/common.hpp"
#include "tbb/parallel_for.h"
#include "tbb/task_arena.h"
#include "util/include/util.hpp"

namespace leonova_a_radix_merge_sort {

tbb::task_arena &LeonovaARadixMergeSortTBB::GetTbbArena() {
  static int num_threads = std::max(1, ppc::util::GetNumThreads());
  static tbb::task_arena arena(num_threads);
  return arena;
}

inline uint64_t LeonovaARadixMergeSortTBB::ToUnsignedValue(int64_t value) {
  return static_cast<uint64_t>(value) ^ kSignBitMask;
}

inline std::pair<size_t, size_t> LeonovaARadixMergeSortTBB::GetChunk(size_t tid, size_t num_threads, size_t size) {
  const size_t chunk = (size + num_threads - 1) / num_threads;
  const size_t begin = tid * chunk;
  const size_t end = std::min(begin + chunk, size);
  return {begin, end};
}

void LeonovaARadixMergeSortTBB::FillUnsignedKeys(const std::vector<int64_t> &arr, size_t left, size_t size,
                                                 std::vector<uint64_t> &keys, size_t num_threads) {
  tbb::parallel_for(static_cast<size_t>(0), num_threads, [&](size_t tid) {
    auto [begin, end] = GetChunk(tid, num_threads, size);
    for (size_t index = begin; index < end; ++index) {
      keys[index] = ToUnsignedValue(arr[left + index]);
    }
  });
}

void LeonovaARadixMergeSortTBB::CountBytesParallel(const std::vector<uint64_t> &keys, size_t size, int shift,
                                                   std::vector<CounterRow> &local_counts, size_t num_threads) {
  tbb::parallel_for(static_cast<size_t>(0), num_threads, [&](size_t tid) {
    auto [begin, end] = GetChunk(tid, num_threads, size);

    auto &row = local_counts[tid];

    for (size_t index = begin; index < end; ++index) {
      ++row[(keys[index] >> shift) & 0xFFU];
    }
  });
}

void LeonovaARadixMergeSortTBB::ReduceCounts(const std::vector<CounterRow> &local_counts, CounterRow &global_counts) {
  std::ranges::fill(global_counts, 0);

  for (const auto &row : local_counts) {
    for (size_t index = 0; index < kNumCounters; ++index) {
      global_counts[index] += row[index];
    }
  }
}

void LeonovaARadixMergeSortTBB::BuildOffsets(const std::vector<CounterRow> &local_counts,
                                             std::vector<CounterRow> &local_offsets, CounterRow &global_counts,
                                             size_t num_threads) {
  CounterRow bucket_totals = global_counts;

  size_t prefix = 0;
  for (size_t index = 0; index < kNumCounters; ++index) {
    size_t count = bucket_totals[index];
    bucket_totals[index] = prefix;
    prefix += count;
  }

  for (size_t tndex = 0; tndex < num_threads; ++tndex) {
    auto &offset_row = local_offsets[tndex];
    const auto &count_row = local_counts[tndex];

    for (size_t index = 0; index < kNumCounters; ++index) {
      offset_row[index] = bucket_totals[index];
      bucket_totals[index] += count_row[index];
    }
  }
}

void LeonovaARadixMergeSortTBB::ScatterParallel(const std::vector<uint64_t> &keys, const std::vector<int64_t> &arr,
                                                size_t left, size_t size, int shift,
                                                std::vector<CounterRow> &local_offsets, std::vector<int64_t> &temp_arr,
                                                std::vector<uint64_t> &temp_keys, size_t num_threads) {
  tbb::parallel_for(static_cast<size_t>(0), num_threads, [&](size_t tid) {
    auto [begin, end] = GetChunk(tid, num_threads, size);

    auto &offsets = local_offsets[tid];

    for (size_t index = begin; index < end; ++index) {
      const size_t byte = (keys[index] >> shift) & 0xFFU;
      const size_t pos = offsets[byte]++;

      temp_arr[pos] = arr[left + index];
      temp_keys[pos] = keys[index];
    }
  });
}

LeonovaARadixMergeSortTBB::LeonovaARadixMergeSortTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<int64_t>(GetInput().size());
}

bool LeonovaARadixMergeSortTBB::ValidationImpl() {
  return !GetInput().empty();
}

bool LeonovaARadixMergeSortTBB::PreProcessingImpl() {
  return true;
}

bool LeonovaARadixMergeSortTBB::RunImpl() {
  if (!ValidationImpl()) {
    return false;
  }

  GetOutput() = GetInput();

  if (GetOutput().size() > 1) {
    RadixMergeSort(GetOutput(), 0, GetOutput().size());
  }

  return true;
}

bool LeonovaARadixMergeSortTBB::PostProcessingImpl() {
  return true;
}

void LeonovaARadixMergeSortTBB::SequentialRadixSort(std::vector<int64_t> &arr, size_t left, size_t right) {
  const size_t size = right - left;
  if (size <= 1) {
    return;
  }

  std::vector<uint64_t> keys(size);
  std::vector<int64_t> temp_arr(size);
  std::vector<uint64_t> temp_keys(size);

  for (size_t index = 0; index < size; ++index) {
    keys[index] = ToUnsignedValue(arr[left + index]);
  }

  for (int byte_pos = 0; byte_pos < kNumBytes; ++byte_pos) {
    const int shift = byte_pos * kByteSize;

    std::vector<size_t> counts(kNumCounters, 0);
    std::vector<size_t> offsets(kNumCounters, 0);

    for (size_t index = 0; index < size; ++index) {
      ++counts[(keys[index] >> shift) & 0xFFU];
    }

    size_t sum = 0;
    for (size_t index = 0; index < kNumCounters; ++index) {
      offsets[index] = sum;
      sum += counts[index];
    }

    for (size_t index = 0; index < size; ++index) {
      size_t byte = (keys[index] >> shift) & 0xFFU;
      size_t pos = offsets[byte]++;

      temp_arr[pos] = arr[left + index];
      temp_keys[pos] = keys[index];
    }

    std::ranges::copy(temp_arr, arr.begin() + static_cast<std::ptrdiff_t>(left));
    keys.swap(temp_keys);
  }
}

void LeonovaARadixMergeSortTBB::RadixSort(std::vector<int64_t> &arr, size_t left, size_t right) {
  const size_t size = right - left;
  if (size <= 1) {
    return;
  }

  if (size < kMinParallelSize) {
    SequentialRadixSort(arr, left, right);
    return;
  }

  auto &arena = GetTbbArena();

  const size_t num_threads = std::max<size_t>(1, std::min<size_t>(arena.max_concurrency(), size));

  std::vector<uint64_t> keys(size);
  std::vector<uint64_t> temp_keys(size);
  std::vector<int64_t> temp_arr(size);

  std::vector<CounterRow> local_counts(num_threads, CounterRow(kNumCounters, 0));
  std::vector<CounterRow> local_offsets(num_threads, CounterRow(kNumCounters, 0));
  CounterRow global_counts(kNumCounters);

  arena.execute([&] {
    FillUnsignedKeys(arr, left, size, keys, num_threads);

    for (int byte_pos = 0; byte_pos < kNumBytes; ++byte_pos) {
      const int shift = byte_pos * kByteSize;

      for (auto &row : local_counts) {
        std::ranges::fill(row, 0);
      }

      CountBytesParallel(keys, size, shift, local_counts, num_threads);
      ReduceCounts(local_counts, global_counts);
      BuildOffsets(local_counts, local_offsets, global_counts, num_threads);

      ScatterParallel(keys, arr, left, size, shift, local_offsets, temp_arr, temp_keys, num_threads);

      std::ranges::copy(temp_arr, arr.begin() + static_cast<std::ptrdiff_t>(left));
      keys.swap(temp_keys);
    }
  });
}

void LeonovaARadixMergeSortTBB::SimpleMerge(std::vector<int64_t> &arr, size_t left, size_t mid, size_t right) {
  std::vector<int64_t> merged(right - left);

  size_t in = left;
  size_t j = mid;
  size_t k = 0;

  while (in < mid && j < right) {
    merged[k++] = (arr[in] <= arr[j]) ? arr[in++] : arr[j++];
  }
  while (in < mid) {
    merged[k++] = arr[in++];
  }
  while (j < right) {
    merged[k++] = arr[j++];
  }

  std::ranges::copy(merged, arr.begin() + static_cast<std::ptrdiff_t>(left));
}

void LeonovaARadixMergeSortTBB::RadixMergeSort(std::vector<int64_t> &arr, size_t left, size_t right) {
  struct Task {
    size_t l, r;
    bool sorted;
  };

  std::vector<Task> stack;
  stack.push_back({left, right, false});

  while (!stack.empty()) {
    auto [l, r, sorted] = stack.back();
    stack.pop_back();

    size_t size = r - l;
    if (size <= 1) {
      continue;
    }

    if (size <= kRadixThreshold) {
      RadixSort(arr, l, r);
      continue;
    }

    size_t mid = l + (size / 2);

    if (!sorted) {
      stack.push_back({l, r, true});
      stack.push_back({mid, r, false});
      stack.push_back({l, mid, false});
    } else {
      SimpleMerge(arr, l, mid, r);
    }
  }
}

}  // namespace leonova_a_radix_merge_sort
