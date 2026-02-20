#include "dorofeev_i_bitwise_sort_double_eo_batcher_merge/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstring>
#include <limits>
#include <vector>

namespace dorofeev_i_bitwise_sort_double_eo_batcher_merge {

namespace {

// Преобразование double в uint64_t с сохранением порядка сортировки
uint64_t DoubleToUint(double d) {
  uint64_t u;
  std::memcpy(&u, &d, sizeof(double));
  if ((u & 0x8000000000000000ULL) != 0) {
    u = ~u; // Инвертируем все биты для отрицательных
  } else {
    u |= 0x8000000000000000ULL; // Инвертируем только знаковый бит для положительных
  }
  return u;
}

// Обратное преобразование
double UintToDouble(uint64_t u) {
  if ((u & 0x8000000000000000ULL) != 0) {
    u &= ~0x8000000000000000ULL;
  } else {
    u = ~u;
  }
  double d;
  std::memcpy(&d, &u, sizeof(double));
  return d;
}

// Поразрядная сортировка подмассива
void RadixSortDouble(std::vector<double>& arr) {
  if (arr.empty()) return;
  
  std::vector<uint64_t> uarr(arr.size());
  for (size_t i = 0; i < arr.size(); ++i) {
    uarr[i] = DoubleToUint(arr[i]);
  }

  std::vector<uint64_t> temp(uarr.size());
  for (int byte = 0; byte < 8; ++byte) {
    int count[256] = {0};
    for (uint64_t val : uarr) {
      count[(val >> (byte * 8)) & 0xFF]++;
    }
    for (int i = 1; i < 256; ++i) {
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

// Четно-нечетное слияние Бэтчера
void OddEvenMerge(std::vector<double>& arr, int l, int n, int step) {
  int m = step * 2;
  if (m < n) {
    OddEvenMerge(arr, l, n, m);
    OddEvenMerge(arr, l + step, n, m);
    for (int i = l + step; i + step < l + n; i += m) {
      if (arr[i] > arr[i + step]) {
        std::swap(arr[i], arr[i + step]);
      }
    }
  } else {
    if (l + step < l + n && arr[l] > arr[l + step]) {
      std::swap(arr[l], arr[l + step]);
    }
  }
}

}  // namespace

DorofeevIBitwiseSortDoubleEOBatcherMergeSEQ::DorofeevIBitwiseSortDoubleEOBatcherMergeSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool DorofeevIBitwiseSortDoubleEOBatcherMergeSEQ::ValidationImpl() {
  return true; // Сортировка пустого массива тоже валидна
}

bool DorofeevIBitwiseSortDoubleEOBatcherMergeSEQ::PreProcessingImpl() {
  local_data_ = GetInput();
  return true;
}

bool DorofeevIBitwiseSortDoubleEOBatcherMergeSEQ::RunImpl() {
  if (local_data_.empty()) return true;

  size_t original_size = local_data_.size();
  
  // Добиваем до ближайшей степени двойки
  size_t pow2 = 1;
  while (pow2 < original_size) pow2 *= 2;
  
  if (pow2 > original_size) {
    local_data_.resize(pow2, std::numeric_limits<double>::max());
  }

  // Для SEQ версии бьем на 2 части, сортируем их и сливаем Бэтчером (имитация параллелизма)
  size_t mid = pow2 / 2;
  std::vector<double> left(local_data_.begin(), local_data_.begin() + mid);
  std::vector<double> right(local_data_.begin() + mid, local_data_.end());

  RadixSortDouble(left);
  RadixSortDouble(right);

  std::copy(left.begin(), left.end(), local_data_.begin());
  std::copy(right.begin(), right.end(), local_data_.begin() + mid);

  OddEvenMerge(local_data_, 0, pow2, 1);

  // Отрезаем фиктивные элементы
  if (pow2 > original_size) {
    local_data_.resize(original_size);
  }

  return true;
}

bool DorofeevIBitwiseSortDoubleEOBatcherMergeSEQ::PostProcessingImpl() {
  GetOutput() = local_data_;
  return true;
}

}  // namespace dorofeev_i_bitwise_sort_double_eo_batcher_merge