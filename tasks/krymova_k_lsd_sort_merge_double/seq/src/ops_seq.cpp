#include "krymova_k_lsd_sort_merge_double/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <vector>

namespace krymova_k_lsd_sort_merge_double {

KrymovaKLsdSortMergeDoubleSEQ::KrymovaKLsdSortMergeDoubleSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType();
}

bool KrymovaKLsdSortMergeDoubleSEQ::ValidationImpl() {
  return true;
}

bool KrymovaKLsdSortMergeDoubleSEQ::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

uint64_t KrymovaKLsdSortMergeDoubleSEQ::DoubleToULL(double d) {
  uint64_t ull = 0;
  std::memcpy(&ull, &d, sizeof(double));

  if ((ull & 0x8000000000000000ULL) != 0U) {
    ull = ~ull;
  } else {
    ull |= 0x8000000000000000ULL;
  }

  return ull;
}

double KrymovaKLsdSortMergeDoubleSEQ::ULLToDouble(uint64_t ull) {
  if ((ull & 0x8000000000000000ULL) != 0U) {
    ull &= 0x7FFFFFFFFFFFFFFFULL;
  } else {
    ull = ~ull;
  }

  double d = 0.0;
  std::memcpy(&d, &ull, sizeof(double));
  return d;
}

void KrymovaKLsdSortMergeDoubleSEQ::LSDSortDouble(double* arr, int size) {
  if (size <= 1) {
    return;
  }

  const int k_bits_per_pass = 8;
  const int k_radix = 1 << k_bits_per_pass;
  const int k_passes = static_cast<int>(sizeof(double)) * 8 / k_bits_per_pass;

  std::vector<uint64_t> ull_arr(size);
  std::vector<uint64_t> ull_tmp(size);

  for (int i = 0; i < size; ++i) {
    ull_arr[i] = DoubleToULL(arr[i]);
  }

  std::vector<unsigned int> count(k_radix, 0U);

  for (int pass = 0; pass < k_passes; ++pass) {
    int shift = pass * k_bits_per_pass;

    std::ranges::fill(count, 0U);

    for (int i = 0; i < size; ++i) {
      unsigned int digit = (ull_arr[i] >> shift) & (k_radix - 1);
      ++count[digit];
    }

    for (int i = 1; i < k_radix; ++i) {
      count[i] += count[i - 1];
    }

    for (int i = size - 1; i >= 0; --i) {
      unsigned int digit = (ull_arr[i] >> shift) & (k_radix - 1);
      ull_tmp[--count[digit]] = ull_arr[i];
    }

    ull_arr.swap(ull_tmp);
  }

  for (int i = 0; i < size; ++i) {
    arr[i] = ULLToDouble(ull_arr[i]);
  }
}

void KrymovaKLsdSortMergeDoubleSEQ::IterativeMergeSort(double* arr, double* tmp, int size, int portion) {
  if (size <= 1) {
    return;
  }

  // ШАГ 1: Сортируем все кусочки размера portion с помощью LSD
  for (int i = 0; i < size; i += portion) {
    int current_size = std::min(portion, size - i);
    LSDSortDouble(arr + i, current_size);
  }

  // ШАГ 2: Итеративно сливаем кусочки
  for (int merge_size = portion; merge_size < size; merge_size *= 2) {
    for (int i = 0; i < size; i += 2 * merge_size) {
      int left_size = merge_size;
      int right_size = std::min(merge_size, size - (i + merge_size));

      if (right_size <= 0) {
        continue;
      }

      double* left = arr + i;
      double* right = arr + i + left_size;
      double* temp = tmp + i;

      std::copy(left, left + left_size, temp);

      int l = 0;
      int r = 0;
      int k = 0;

      while (l < left_size && r < right_size) {
        if (temp[l] <= right[r]) {
          left[k++] = temp[l++];
        } else {
          left[k++] = right[r++];
        }
      }

      while (l < left_size) {
        left[k++] = temp[l++];
      }
    }
  }
}

bool KrymovaKLsdSortMergeDoubleSEQ::RunImpl() {
  OutType& output = GetOutput();
  int size = static_cast<int>(output.size());

  if (size <= 1) {
    return true;
  }

  std::vector<double> tmp(size);
  int portion = std::max(1, size / 10);

  IterativeMergeSort(output.data(), tmp.data(), size, portion);

  return true;
}

bool KrymovaKLsdSortMergeDoubleSEQ::PostProcessingImpl() {
  const OutType& output = GetOutput();

  for (size_t i = 1; i < output.size(); ++i) {
    if (output[i] < output[i - 1]) {
      return false;
    }
  }

  return true;
}

}  // namespace krymova_k_lsd_sort_merge_double