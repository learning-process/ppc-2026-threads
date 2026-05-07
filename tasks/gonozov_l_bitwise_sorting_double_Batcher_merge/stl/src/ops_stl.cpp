#include "gonozov_l_bitwise_sorting_double_Batcher_merge/STL/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <execution>
#include <future>
#include <limits>
#include <vector>

#include "gonozov_l_bitwise_sorting_double_Batcher_merge/common/include/common.hpp"

namespace gonozov_l_bitwise_sorting_double_batcher_merge {

GonozovLBitSortBatcherMergeSTL::GonozovLBitSortBatcherMergeSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool GonozovLBitSortBatcherMergeSTL::ValidationImpl() {
  return !GetInput().empty();  // проверка на то, что исходный массив непустой
}

bool GonozovLBitSortBatcherMergeSTL::PreProcessingImpl() {
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

void RadixSortDouble(std::vector<double>& data,
                     size_t begin,
                     size_t end) {
  if (end <= begin + 1) {
    return;
  }

  size_t size = end - begin;

  std::vector<uint64_t> keys(size);

  for (size_t i = 0; i < size; ++i) {
    keys[i] = DoubleToSortableInt(data[begin + i]);
  }

  constexpr int radix = 256;

  std::vector<uint64_t> temp_keys(size);

  for (int pass = 0; pass < 8; ++pass) {

    std::array<size_t, radix> count{};

    int shift = pass * 8;

    for (uint64_t key : keys) {
      count[(key >> shift) & 0xFF]++;
    }

    for (int i = 1; i < radix; ++i) {
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

void MergingHalves(std::vector<double> &arr, size_t i, size_t len) {  // слияние половинок
  size_t half = len / 2;
  size_t end = std::min(i + len, arr.size());

  for (size_t step = half; step > 0; step /= 2) {
    for (size_t j = i; j + step < end; ++j) {
      if (arr[j] > arr[j + step]) {
        std::swap(arr[j], arr[j + step]);
      }
    }
  }
}

// void BatcherOddEvenMergeIterative(std::vector<double> &arr, size_t n) {
//   if (n <= 1) {
//     return;
//   }
//   n = std::min(n, arr.size());
//   // Сначала сливаем блоки размером 1, потом 2, потом 4 и т.д.
//   for (size_t len = 2; len <= n; len *= 2) {
//     for (size_t i = 0; i < n; i += len) {
//       MergingHalves(arr, i, len);
//     }
//   }
// }
// void BatcherOddEvenMergeIterative(std::vector<double>& arr, size_t n) {
//   if (n <= 1) {
//     return;
//   }

//   n = std::min(n, arr.size());

//   // Слияние блоков: 2,4,8,16...
//   for (size_t len = 2; len <= n; len *= 2) {

//     // Индексы независимых блоков
//     std::vector<size_t> blocks;
//     blocks.reserve((n + len - 1) / len);

//     for (size_t i = 0; i < n; i += len) {
//       blocks.push_back(i);
//     }

//     // Параллельное выполнение merge блоков
//     std::for_each(std::execution::par,
//                   blocks.begin(),
//                   blocks.end(),
//                   [&](size_t i) {
//                     MergingHalves(arr, i, len);
//                   });
//   }
// }

// Нахождение ближайшей степени двойки, большей или равной n
size_t NextPowerOfTwo(size_t n) {
  size_t power = 1;
  while (power < n) {
    power <<= 1;
  }
  return power;
}

void HybridSortDouble(std::vector<double>& data) {
  if (data.size() <= 1) {
    return;
  }

  size_t original_size = data.size();

  size_t new_size = NextPowerOfTwo(original_size);

  data.resize(
      new_size,
      std::numeric_limits<double>::infinity());

  size_t threads =
      std::thread::hardware_concurrency();

  size_t block_size =
      (new_size + threads - 1) / threads;

  // PARALLEL RADIX SORT

  std::vector<std::future<void>> futures;

  futures.reserve(threads);

  for (size_t t = 0; t < threads; ++t) {

    futures.push_back(
        std::async(std::launch::async,
        [&, t]() {

          size_t begin = t * block_size;

          size_t end =
              std::min(begin + block_size,
                       new_size);

          if (begin < end) {
            RadixSortDouble(data, begin, end);
          }
        }));
  }

  for (auto& f : futures) {
    f.get();
  }

  // PARALLEL BATCHER MERGE

  for (size_t merge_size = block_size;
       merge_size < new_size;
       merge_size *= 2) {

    size_t blocks_count =
        (new_size + merge_size * 2 - 1) /
        (merge_size * 2);

    std::vector<size_t> indices(blocks_count);

    std::iota(indices.begin(),
              indices.end(),
              0);

    std::for_each(std::execution::par,
                  indices.begin(),
                  indices.end(),
                  [&](size_t block_id) {

      size_t begin =
          block_id * merge_size * 2;

      size_t len =
          std::min(
              merge_size * 2,
              new_size - begin);

      if (begin < new_size) {
        MergingHalves(
            data,
            begin,
            len);
      }
    });
  }

  data.resize(original_size);
}

}  // namespace

bool GonozovLBitSortBatcherMergeSTL::RunImpl() {
  std::vector<double> array = GetInput();
  HybridSortDouble(array);
  GetOutput() = array;
  return true;
}

bool GonozovLBitSortBatcherMergeSTL::PostProcessingImpl() {
  return true;
}

}  // namespace gonozov_l_bitwise_sorting_double_batcher_merge
