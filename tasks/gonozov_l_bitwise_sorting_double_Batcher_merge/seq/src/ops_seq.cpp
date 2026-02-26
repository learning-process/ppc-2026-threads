#include "gonozov_l_bitwise_sorting_double_Batcher_merge/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>

#include "gonozov_l_bitwise_sorting_double_Batcher_merge/common/include/common.hpp"

namespace gonozov_l_bitwise_sorting_double_batcher_merge {

GonozovLBitSortBatcherMergeSEQ::GonozovLBitSortBatcherMergeSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool GonozovLBitSortBatcherMergeSEQ::ValidationImpl() {
  return !GetInput().empty();  // проверка на то, что исходный массив непустой
}

bool GonozovLBitSortBatcherMergeSEQ::PreProcessingImpl() {
  return true;
}

namespace {
/// double -> uint64_t
uint64_t DoubleToSortableInt(double d) {
  uint64_t bits = 0;
  std::memcpy(&bits, &d, sizeof(double));
  if ((bits >> 63) != 0) {  // Отрицательное число
    return ~bits;           // Инвертируем ВСЕ биты
  }  // Положительное число или ноль
  return bits ^ 0x8000000000000000ULL;  // Инвертируем ТОЛЬКО знаковый бит
}

// uint64_t -> double
double SortableIntToDouble(uint64_t bits) {
  if ((bits >> 63) != 0) {                // Если старший бит установлен (было положительное)
    bits = bits ^ 0x8000000000000000ULL;  // Возвращаем знаковый бит
  } else {                                // Если старший бит не установлен (было отрицательное)
    bits = ~bits;                         // Инвертируем все биты обратно
  }

  double result = NAN;
  std::memcpy(&result, &bits, sizeof(double));
  return result;
}

void RadixSortDouble(std::vector<double> &data) {
  if (data.empty()) {
    return;
  }

  // Преобразуем в сортируемые целые числа
  std::vector<uint64_t> keys(data.size());
  for (size_t i = 0; i < data.size(); ++i) {
    keys[i] = DoubleToSortableInt(data[i]);
  }

  const int radix = 256;  // 8 бит за проход
  std::vector<uint64_t> temp_keys(data.size());

  // 8 проходов для 64-битных чисел (8 байт)
  for (int pass = 0; pass < 8; ++pass) {
    std::vector<size_t> count(radix, 0);
    int shift = pass * 8;
    // Подсчет
    for (uint64_t key : keys) {
      uint8_t byte = (key >> shift) & 0xFF;
      count[byte]++;
    }

    // Накопление
    for (int i = 1; i < radix; ++i) {
      count[i] += count[i - 1];
    }

    // Распределение
    for (int i = static_cast<int>(keys.size() - 1); i >= 0; --i) {
      uint8_t byte = (keys[i] >> shift) & 0xFF;
      temp_keys[--count[byte]] = keys[i];
    }
    keys.swap(temp_keys);
  }

  // Преобразуем обратно
  for (size_t i = 0; i < data.size(); ++i) {
    data[i] = SortableIntToDouble(keys[i]);
  }
}

void BatcherOddEvenMergeIterative(std::vector<double> &arr, size_t n) {
  if (n <= 1) {
    return;
  }

  for (size_t block_size = 1; block_size < n; block_size *= 2) {
    for (size_t step = block_size; step > 0; step /= 2) {
      for (size_t i = 0; i < n - step; ++i) {
        // Проверяем, нужно ли сравнивать элементы i и i+step
        if ((i & block_size) == 0) {
          size_t idx1 = i;
          size_t idx2 = i + step;
          if (idx2 < n && arr[idx1] > arr[idx2]) {
            std::swap(arr[idx1], arr[idx2]);
          }
        }
      }
    }
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
  data.resize(new_size, std::numeric_limits<double>::max());

  size_t mid = new_size / 2;
  std::vector<double> left(data.begin(), data.begin() + mid);
  std::vector<double> right(data.begin() + mid, data.end());

  // Сортируем каждую половину поразрядно
  RadixSortDouble(left);
  RadixSortDouble(right);

  // Собираем обратно в единый массив
  std::copy(left.begin(), left.end(), data.begin());
  std::copy(right.begin(), right.end(), data.begin() + mid);

  // Используем слияние Бэтчера для слияния двух отсортированных массивов
  BatcherOddEvenMergeIterative(data, new_size);

  // Обрезаем до исходного размера
  data.resize(original_size);
}

}  // namespace

bool GonozovLBitSortBatcherMergeSEQ::RunImpl() {
  std::vector<double> array = GetInput();
  HybridSortDouble(array);
  GetOutput() = array;
  return true;
}

bool GonozovLBitSortBatcherMergeSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace gonozov_l_bitwise_sorting_double_batcher_merge
