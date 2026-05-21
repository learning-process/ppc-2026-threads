#include "gonozov_l_bitwise_sorting_double_Batcher_merge/tbb/include/ops_tbb.hpp"

#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/parallel_invoke.h>
#include <tbb/tbb.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

#include "gonozov_l_bitwise_sorting_double_Batcher_merge/common/include/common.hpp"

namespace gonozov_l_bitwise_sorting_double_batcher_merge {

GonozovLBitSortBatcherMergeTBB::GonozovLBitSortBatcherMergeTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool GonozovLBitSortBatcherMergeTBB::ValidationImpl() {
  return !GetInput().empty();  // проверка на то, что исходный массив непустой
}

bool GonozovLBitSortBatcherMergeTBB::PreProcessingImpl() {
  return true;
}

namespace {
/// double -> uint64_t
uint64_t DoubleToSortableInt(double d) {
  uint64_t bits = 0;
  std::memcpy(&bits, &d, sizeof(double));
  if ((bits >> 63) != 0) {  // Отрицательное число
    return ~bits;           // Инвертируем все биты
  }  // Положительное число или ноль
  return bits | 0x8000000000000000ULL;
}

// uint64_t -> double
double SortableIntToDouble(uint64_t bits) {
  if ((bits >> 63) != 0) {           // Если старший бит установлен (было положительное)
    bits &= ~0x8000000000000000ULL;  // Убираем старший бит
  } else {                           // Если старший бит не установлен (было отрицательное число)
    bits = ~bits;                    // Инвертируем все биты обратно
  }

  double result = 0.0;
  std::memcpy(&result, &bits, sizeof(double));
  return result;
}

inline size_t NextPowerOfTwo(size_t n) {
  size_t power = 1;

  while (power < n) {
    power <<= 1;
  }

  return power;
}

constexpr size_t kRadix = 256;

void RadixSortDoubleChunk(std::vector<uint64_t> &data, size_t start, size_t end) {
  size_t size = end - start;
  if (size <= 1) {
    return;
  }

  std::vector<uint64_t> temp(size);
  uint64_t *keys = data.data() + start;
  uint64_t *temp_ptr = temp.data();

  for (int pass = 0; pass < 8; ++pass) {
    std::array<size_t, kRadix> count{};
    int shift = pass * 8;

    for (size_t i = 0; i < size; ++i) {
      ++count[(keys[i] >> shift) & 0xFF];
    }

    size_t sum = 0;
    for (size_t i = 0; i < kRadix; ++i) {
      size_t c = count[i];
      count[i] = sum;
      sum += c;
    }

    for (size_t i = 0; i < size; ++i) {
      uint8_t byte = (keys[i] >> shift) & 0xFF;
      temp_ptr[count[byte]++] = keys[i];
    }

    std::memcpy(keys, temp_ptr, size * sizeof(uint64_t));
  }
}

void CompareExchangeBlocks(uint64_t *arr, size_t i, size_t step) {
  for (size_t k = 0; k < step; ++k) {
    if (arr[i + k] > arr[i + k + step]) {
      std::swap(arr[i + k], arr[i + k + step]);
    }
  }

void OddEvenMergeIterative(uint64_t *arr, size_t start, size_t n) {
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

void ProcessChunkTBB(std::vector<uint64_t> &keys, int chunk_idx, size_t chunk_size) {
  size_t start_idx = static_cast<size_t>(chunk_idx) * chunk_size;
  RadixSortDoubleChunk(keys, start_idx, start_idx + chunk_size);
}

void TBBSort(std::vector<uint64_t> &keys, size_t pow2, size_t chunk_size, int num_chunks_int, int num_threads) {
  tbb::task_arena arena(num_threads);
  arena.execute([&]() {
    tbb::parallel_for(tbb::blocked_range<int>(0, num_chunks_int), [&](const tbb::blocked_range<int> &r) {
      for (int i = r.begin(); i != r.end(); ++i) {
        ProcessChunkTBB(keys, i, chunk_size);
      }
    });

    for (size_t size = chunk_size; size < pow2; size *= 2) {
      int merges_count = static_cast<int>(pow2 / (size * 2));

      tbb::parallel_for(tbb::blocked_range<int>(0, merges_count), [&](const tbb::blocked_range<int> &r) {
        for (int i = r.begin(); i != r.end(); ++i) {
          OddEvenMergeIterative(keys.data(), static_cast<size_t>(i) * 2 * size, 2 * size);
        }
      });
    }
  });
}

void HybridSort(std::vector<double> &arr) {
  if (arr.empty()) {
    return;
  }

  size_t original_size = arr.size();
  size_t pow2 = NextPowerOfTwo(original_size);

  std::vector<uint64_t> keys(pow2);
  for (size_t i = 0; i < original_size; ++i) {
    keys[i] = DoubleToSortableInt(arr[i]);
  }

  uint64_t max_key = DoubleToSortableInt(std::numeric_limits<double>::max());
  for (size_t i = original_size; i < pow2; ++i) {
    keys[i] = max_key;
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
  int num_chunks_int = static_cast<int>(num_chunks);

  TBBSort(keys, pow2, chunk_size, num_chunks_int, num_threads);

  for (size_t i = 0; i < original_size; ++i) {
    arr[i] = SortableIntToDouble(keys[i]);
  }
}

}  // namespace

bool GonozovLBitSortBatcherMergeTBB::RunImpl() {
  std::vector<double> array = GetInput();
  HybridSort(array);
  GetOutput() = array;
  return true;
}

bool GonozovLBitSortBatcherMergeTBB::PostProcessingImpl() {
  return true;
}

}  // namespace gonozov_l_bitwise_sorting_double_batcher_merge
