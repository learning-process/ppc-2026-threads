#include "gonozov_l_bitwise_sorting_double_Batcher_merge/all/include/ops_all.hpp"

#include <mpi.h>
#include <tbb/tbb.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

#include "gonozov_l_bitwise_sorting_double_Batcher_merge/common/include/common.hpp"
#include "util/include/util.hpp"

namespace gonozov_l_bitwise_sorting_double_batcher_merge {

namespace {

uint64_t DoubleToSortableInt(double d) {
  uint64_t bits = 0;
  std::memcpy(&bits, &d, sizeof(double));

  if ((bits >> 63) != 0) {
    return ~bits;
  }

  return bits | 0x8000000000000000ULL;
}

double SortableIntToDouble(uint64_t bits) {
  if ((bits >> 63) != 0) {
    bits &= ~0x8000000000000000ULL;
  } else {
    bits = ~bits;
  }

  double result = 0.0;
  std::memcpy(&result, &bits, sizeof(double));

  return result;
}

size_t NextPowerOfTwo(size_t n) {
  size_t power = 1;
  while (power < n) {
    power <<= 1;
  }
  return power;
}

void RadixSortDouble(std::vector<double> &data) {
  if (data.empty()) {
    return;
  }

  std::vector<uint64_t> keys(data.size());

  for (size_t i = 0; i < data.size(); ++i) {
    keys[i] = DoubleToSortableInt(data[i]);
  }

  std::vector<uint64_t> temp(keys.size());

  constexpr int kRadix = 256;

  for (int pass = 0; pass < 8; ++pass) {
    std::vector<int> count(kRadix, 0);

    int shift = pass * 8;

    for (uint64_t key : keys) {
      count[(key >> shift) & 0xFF]++;
    }

    for (int i = 1; i < kRadix; ++i) {
      count[i] += count[i - 1];
    }

    for (int i = static_cast<int>(keys.size()) - 1; i >= 0; --i) {
      uint8_t byte = (keys[i] >> shift) & 0xFF;
      temp[--count[byte]] = keys[i];
    }

    std::swap(keys, temp);
  }

  for (size_t i = 0; i < data.size(); ++i) {
    data[i] = SortableIntToDouble(keys[i]);
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

void SortChunkALL(double *raw_data, int chunk_idx, size_t chunk_size) {
  size_t start_idx = static_cast<size_t>(chunk_idx) * chunk_size;

  std::vector<double> local_arr(raw_data + start_idx, raw_data + start_idx + chunk_size);

  RadixSortDouble(local_arr);

  std::ranges::copy(local_arr.begin(), local_arr.end(), raw_data + start_idx);
}

}  // namespace

GonozovLBitSortBatcherMergeALL::GonozovLBitSortBatcherMergeALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());

  GetInput() = in;
}

bool GonozovLBitSortBatcherMergeALL::ValidationImpl() {
  return true;
}

bool GonozovLBitSortBatcherMergeALL::PreProcessingImpl() {
  local_data_ = GetInput();

  return true;
}

bool GonozovLBitSortBatcherMergeALL::RunImpl() {
  if (local_data_.empty()) {
    return true;
  }

  size_t n = local_data_.size();

  size_t new_size = NextPowerOfTwo(n);

  if (new_size > n) {
    local_data_.resize(new_size, std::numeric_limits<double>::infinity());
  }

  int num_threads = ppc::util::GetNumThreads();

  if (num_threads <= 0) {
    num_threads = 1;
  }

  size_t num_chunks = 1;

  while (num_chunks * 2 <= static_cast<size_t>(num_threads) && num_chunks * 2 <= new_size) {
    num_chunks *= 2;
  }

  size_t chunk_size = new_size / num_chunks;
  chunk_size = std::max<size_t>(1, chunk_size);
  double *raw_data = local_data_.data();
  int num_chunks_int = static_cast<int>(num_chunks);

  tbb::parallel_for(0, num_chunks_int, [&](int i) { SortChunkALL(raw_data, i, chunk_size); });

  for (size_t size = chunk_size; size < new_size; size *= 2) {
    int merges_count = static_cast<int>(new_size / (size * 2 + 0.0000000001));

    tbb::parallel_for(0, merges_count,
                      [&](int i) { OddEvenMergeIterative(raw_data, static_cast<size_t>(i) * 2 * size, 2 * size); });
  }

  if (new_size > n) {
    local_data_.resize(n);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  return true;
}

bool GonozovLBitSortBatcherMergeALL::PostProcessingImpl() {
  GetOutput() = local_data_;

  return true;
}

}  // namespace gonozov_l_bitwise_sorting_double_batcher_merge
