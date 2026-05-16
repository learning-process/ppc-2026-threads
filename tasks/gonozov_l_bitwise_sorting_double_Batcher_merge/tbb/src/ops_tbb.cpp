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

void RadixSortDouble(std::vector<double> &data, size_t begin, size_t end) {
  if (end <= begin + 1) {
    return;
  }

  size_t size = end - begin;

  std::vector<uint64_t> keys(size);

  for (size_t i = 0; i < size; ++i) {
    keys[i] = DoubleToSortableInt(data[begin + i]);
  }

  constexpr int kRadix = 256;

  std::vector<uint64_t> temp_keys(size);

  for (int pass = 0; pass < 8; ++pass) {
    std::vector<size_t> count(kRadix, 0);

    int shift = pass * 8;

    for (uint64_t key : keys) {
      count[(key >> shift) & 0xFF]++;
    }

    for (int i = 1; i < kRadix; ++i) {
      count[i] += count[i - 1];
    }

    for (int i = static_cast<int>(size) - 1; i >= 0; --i) {
      uint8_t byte = (keys[i] >> shift) & 0xFF;
      temp_keys[--count[byte]] = keys[i];
    }

    keys.swap(temp_keys);
  }

  for (size_t i = 0; i < size; ++i) {
    data[begin + i] = SortableIntToDouble(keys[i]);
  }
}

// Нахождение ближайшей степени двойки, большей или равной n
size_t NextPowerOfTwo(size_t n) {
  size_t power = 1;
  while (power < n) {
    power <<= 1;
  }
  return power;
}

void HybridSortDouble(std::vector<double> &data) {
  if (data.size() <= 1) {
    return;
  }

  size_t original_size = data.size();

  size_t new_size = NextPowerOfTwo(original_size);

  data.resize(new_size, std::numeric_limits<double>::infinity());

  size_t threads = oneapi::tbb::this_task_arena::max_concurrency();

  size_t block_size = NextPowerOfTwo((new_size + threads - 1) / threads);

  // parallel radix sort
  oneapi::tbb::parallel_for(static_cast<size_t>(0), threads, [&](size_t t) {
    size_t begin = t * block_size;

    size_t end = std::min(begin + block_size, new_size);

    if (begin < end) {
      RadixSortDouble(data, begin, end);
    }
  });

  // parallel batcher merge tree
  for (size_t merge_size = block_size; merge_size < new_size; merge_size *= 2) {
    oneapi::tbb::parallel_for(static_cast<size_t>(0), new_size, merge_size * 2, [&](size_t i) {
      size_t mid = std::min(i + merge_size, new_size);
      size_t end = std::min(i + (merge_size * 2), new_size);

      using DiffT = std::vector<double>::difference_type;

      std::inplace_merge(data.begin() + static_cast<DiffT>(i), data.begin() + static_cast<DiffT>(mid),
                         data.begin() + static_cast<DiffT>(end));
    });
  }

  data.resize(original_size);
}

}  // namespace

bool GonozovLBitSortBatcherMergeTBB::RunImpl() {
  std::vector<double> array = GetInput();
  HybridSortDouble(array);
  GetOutput() = array;
  return true;
}

bool GonozovLBitSortBatcherMergeTBB::PostProcessingImpl() {
  return true;
}

}  // namespace gonozov_l_bitwise_sorting_double_batcher_merge
