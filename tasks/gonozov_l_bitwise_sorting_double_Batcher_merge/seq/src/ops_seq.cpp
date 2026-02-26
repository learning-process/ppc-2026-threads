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
    for (size_t i = keys.size() - 1; i >= 0; --i) {
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

// NOLINTNEXTLINE(misc-no-recursion)
void BatcherOddEvenMerge(std::vector<double> &arr, int low, int high) {
  if (high - low <= 1) {
    return;
  }

  int mid = (low + high) / 2;
  BatcherOddEvenMerge(arr, low, mid);
  BatcherOddEvenMerge(arr, mid, high);

  // Сравниваем и меняем элементы из двух половин
  for (int i = low; i < mid; ++i) {
    int j = i + mid - low;  // Индекс во второй половине
    if (j < high && arr[i] > arr[j]) {
      std::swap(arr[i], arr[j]);
    }
  }
}

void HybridSortDouble(std::vector<double> &data) {
  RadixSortDouble(data);
  BatcherOddEvenMerge(data, 0, static_cast<int>(data.size()));
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
