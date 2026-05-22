#include "gonozov_l_bitwise_sorting_double_Batcher_merge/tbb/include/ops_tbb.hpp"

#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/parallel_invoke.h>
#include <tbb/tbb.h>

#include <algorithm>
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
  return !GetInput().empty();
}

bool GonozovLBitSortBatcherMergeTBB::PreProcessingImpl() {
  return true;
}

namespace {

/// double -> uint64_t
uint64_t DoubleToSortableInt(double d) {
  uint64_t bits = 0;
  std::memcpy(&bits, &d, sizeof(double));
  if ((bits >> 63) != 0) {
    return ~bits;
  }
  return bits | 0x8000000000000000ULL;
}

/// uint64_t -> double
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

constexpr size_t kRadix = 256;

/// Поразрядная сортировка для vector<double> (работает с копией)
void RadixSortDouble(std::vector<double> &data) {
  if (data.empty()) {
    return;
  }

  std::vector<uint64_t> keys(data.size());
  for (size_t i = 0; i < data.size(); ++i) {
    keys[i] = DoubleToSortableInt(data[i]);
  }

  std::vector<uint64_t> temp_keys(data.size());

  for (int pass = 0; pass < 8; ++pass) {
    std::vector<size_t> count(kRadix, 0);
    int shift = pass * 8;

    for (uint64_t key : keys) {
      uint8_t byte = (key >> shift) & 0xFF;
      count[byte]++;
    }

    for (int i = 1; std::cmp_less(i, kRadix); ++i) {
      count[i] += count[i - 1];
    }

    for (int i = static_cast<int>(keys.size()) - 1; i >= 0; --i) {
      uint8_t byte = (keys[i] >> shift) & 0xFF;
      temp_keys[--count[byte]] = keys[i];
    }
    keys.swap(temp_keys);
  }

  for (size_t i = 0; i < data.size(); ++i) {
    data[i] = SortableIntToDouble(keys[i]);
  }
}

/// Слияние половинок для сети Бэтчера
void MergingHalves(std::vector<double> &arr, size_t i, size_t len) {
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

/// Итеративная сеть Бэтчера (параллельная)
void BatcherOddEvenMergeIterative(std::vector<double> &arr, size_t n) {
  if (n <= 1) {
    return;
  }
  n = std::min(n, arr.size());

  for (size_t len = 2; len <= n; len *= 2) {
    tbb::parallel_for(tbb::blocked_range<size_t>(0, n, len), [&](const tbb::blocked_range<size_t> &r) {
      for (size_t i = r.begin(); i < r.end(); i += len) {
        MergingHalves(arr, i, len);
      }
    });
  }
}

/// Гибридная сортировка (с копированием половин, как в OMP)
void HybridSortDouble(std::vector<double> &data) {
  if (data.size() <= 1) {
    return;
  }

  size_t original_size = data.size();
  size_t new_size = NextPowerOfTwo(original_size);
  data.resize(new_size, std::numeric_limits<double>::max());

  size_t mid = new_size / 2;

  // Копируем половины
  std::vector<double> left(data.begin(), data.begin() + static_cast<ptrdiff_t>(mid));
  std::vector<double> right(data.begin() + static_cast<ptrdiff_t>(mid), data.end());

  // Параллельно сортируем две половины
  tbb::parallel_invoke([&]() { RadixSortDouble(left); }, [&]() { RadixSortDouble(right); });

  // Собираем обратно
  std::ranges::copy(left.begin(), left.end(), data.begin());
  std::ranges::copy(right.begin(), right.end(), data.begin() + static_cast<ptrdiff_t>(mid));

  // Слияние Бэтчера
  BatcherOddEvenMergeIterative(data, new_size);

  // Обрезаем до исходного размера
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
