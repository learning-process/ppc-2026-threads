#include "dorofeev_i_bitwise_sort_double_eo_batcher_merge/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits>
#include <thread>
#include <utility>
#include <vector>

#include "dorofeev_i_bitwise_sort_double_eo_batcher_merge/common/include/common.hpp"
#include "util/include/util.hpp"

namespace dorofeev_i_bitwise_sort_double_eo_batcher_merge {

namespace {

uint64_t DoubleToUint(double d) {
  uint64_t u = 0;
  std::memcpy(&u, &d, sizeof(double));
  if ((u & 0x8000000000000000ULL) != 0) {
    u = ~u;
  } else {
    u |= 0x8000000000000000ULL;
  }
  return u;
}

double UintToDouble(uint64_t u) {
  if ((u & 0x8000000000000000ULL) != 0) {
    u &= ~0x8000000000000000ULL;
  } else {
    u = ~u;
  }
  double d = 0.0;
  std::memcpy(&d, &u, sizeof(double));
  return d;
}

void RadixSortDouble(std::vector<double> &arr) {
  if (arr.empty()) {
    return;
  }

  std::vector<uint64_t> uarr(arr.size());
  for (size_t i = 0; i < arr.size(); ++i) {
    uarr[i] = DoubleToUint(arr[i]);
  }

  std::vector<uint64_t> temp(uarr.size());
  for (size_t byte = 0; byte < 8; ++byte) {
    std::vector<int> count(256, 0);
    for (uint64_t val : uarr) {
      count[(val >> (byte * 8)) & 0xFF]++;
    }
    for (size_t i = 1; i < 256; ++i) {
      count[i] += count[i - 1];
    }
    for (int i = static_cast<int>(uarr.size()) - 1; i >= 0; --i) {
      temp[--count[(uarr[i] >> (byte * 8)) & 0xFF]] = uarr[i];
    }
    uarr = temp;
  }

  for (size_t i = 0; i < arr.size(); ++i) {
    arr[i] = UintToDouble(uarr[i]);
  }
}

void CompareExchangeBlocks(double *arr, size_t i, size_t step) {
  for (size_t k = 0; k < step; ++k) {
    if (arr[i + k] > arr[i + k + step]) {
      std::swap(arr[i + k], arr[i + k + step]);
    }
  }
}

void OddEvenMergeIterative(double *arr, size_t start, size_t n) {
  if (n <= 1) {
    return;
  }

  size_t step = n / 2;
  CompareExchangeBlocks(arr, start, step);

  step /= 2;
  for (; step > 0; step /= 2) {
    for (size_t i = step; i < n - step; i += step * 2) {
      CompareExchangeBlocks(arr, start + i, step);
    }
  }
}

void ProcessChunkSTL(double *raw_data, int chunk_idx, size_t chunk_size) {
  size_t start_idx = static_cast<size_t>(chunk_idx) * chunk_size;

  std::vector<double> local_arr(chunk_size);
  double *local_raw = local_arr.data();

  for (size_t j = 0; j < chunk_size; ++j) {
    local_raw[j] = raw_data[start_idx + j];
  }

  RadixSortDouble(local_arr);

  for (size_t j = 0; j < chunk_size; ++j) {
    raw_data[start_idx + j] = local_raw[j];
  }
}

void ExecuteSTLSort(double *raw_data, size_t pow2, size_t chunk_size, int num_chunks_int) {
  std::vector<std::thread> threads;
  threads.reserve(num_chunks_int);

  // 1. Параллельно сортируем каждый блок, раздавая задачи сырым std::thread
  for (int i = 0; i < num_chunks_int; ++i) {
    threads.emplace_back([raw_data, i, chunk_size]() { ProcessChunkSTL(raw_data, i, chunk_size); });
  }

  for (auto &t : threads) {
    t.join();
  }
  threads.clear();

  // 2. Собираем отсортированные блоки вместе (Параллельное дерево слияния Бэтчера)
  for (size_t size = chunk_size; size < pow2; size *= 2) {
    int merges_count = static_cast<int>(pow2 / (size * 2));
    threads.reserve(merges_count);

    for (int i = 0; i < merges_count; ++i) {
      threads.emplace_back(
          [raw_data, i, size]() { OddEvenMergeIterative(raw_data, static_cast<size_t>(i) * 2 * size, 2 * size); });
    }

    for (auto &t : threads) {
      t.join();
    }
    threads.clear();
  }
}

}  // namespace

DorofeevIBitwiseSortDoubleEOBatcherMergeSTL::DorofeevIBitwiseSortDoubleEOBatcherMergeSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool DorofeevIBitwiseSortDoubleEOBatcherMergeSTL::ValidationImpl() {
  return true;
}

bool DorofeevIBitwiseSortDoubleEOBatcherMergeSTL::PreProcessingImpl() {
  local_data_ = GetInput();
  return true;
}

bool DorofeevIBitwiseSortDoubleEOBatcherMergeSTL::RunImpl() {
  if (local_data_.empty()) {
    return true;
  }

  size_t original_size = local_data_.size();
  size_t pow2 = 1;
  while (pow2 < original_size) {
    pow2 *= 2;
  }

  if (pow2 > original_size) {
    local_data_.resize(pow2, std::numeric_limits<double>::max());
  }

  int num_threads = ppc::util::GetNumThreads();
  if (num_threads <= 0) {
    num_threads = 1;
  }

  size_t num_chunks = 1;
  while (num_chunks * 2 <= static_cast<size_t>(num_threads) && num_chunks * 2 <= pow2) {
    num_chunks *= 2;
  }

  size_t chunk_size = pow2 / num_chunks;
  double *raw_data = local_data_.data();
  int num_chunks_int = static_cast<int>(num_chunks);

  ExecuteSTLSort(raw_data, pow2, chunk_size, num_chunks_int);

  if (pow2 > original_size) {
    local_data_.resize(original_size);
  }

  return true;
}

bool DorofeevIBitwiseSortDoubleEOBatcherMergeSTL::PostProcessingImpl() {
  GetOutput() = local_data_;
  return true;
}

}  // namespace dorofeev_i_bitwise_sort_double_eo_batcher_merge
