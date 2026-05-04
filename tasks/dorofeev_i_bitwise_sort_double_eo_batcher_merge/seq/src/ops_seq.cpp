#include "dorofeev_i_bitwise_sort_double_eo_batcher_merge/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

#include "dorofeev_i_bitwise_sort_double_eo_batcher_merge/common/include/common.hpp"

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

// Вынесли операцию сравнения и перестановки блоков в отдельную функцию
void CompareExchangeBlocks(std::vector<double> &arr, size_t i, size_t step) {
  for (size_t k = 0; k < step; ++k) {
    if (arr[i + k] > arr[i + k + step]) {
      std::swap(arr[i + k], arr[i + k + step]);
    }
  }
}

// Теперь функция выглядит максимально просто и элегантно
void OddEvenMergeIterative(std::vector<double> &arr, size_t n) {
  if (n <= 1) {
    return;
  }

  // Первая фаза слияния
  size_t step = n / 2;
  CompareExchangeBlocks(arr, 0, step);

  // Последующие фазы
  step /= 2;
  for (; step > 0; step /= 2) {
    for (size_t i = step; i < n - step; i += step * 2) {
      CompareExchangeBlocks(arr, i, step);
    }
  }
}

}  // namespace

DorofeevIBitwiseSortDoubleEOBatcherMergeSEQ::DorofeevIBitwiseSortDoubleEOBatcherMergeSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool DorofeevIBitwiseSortDoubleEOBatcherMergeSEQ::ValidationImpl() {
  return true;
}

bool DorofeevIBitwiseSortDoubleEOBatcherMergeSEQ::PreProcessingImpl() {
  local_data_ = GetInput();
  return true;
}

bool DorofeevIBitwiseSortDoubleEOBatcherMergeSEQ::RunImpl() {
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

  size_t mid = pow2 / 2;
  std::vector<double> left(local_data_.begin(), local_data_.begin() + static_cast<ptrdiff_t>(mid));
  std::vector<double> right(local_data_.begin() + static_cast<ptrdiff_t>(mid), local_data_.end());

  RadixSortDouble(left);
  RadixSortDouble(right);

  std::ranges::copy(left, local_data_.begin());
  std::ranges::copy(right, local_data_.begin() + static_cast<ptrdiff_t>(mid));

  OddEvenMergeIterative(local_data_, pow2);

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
