#include "khruev_a_radix_sorting_int_bather_merge/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "khruev_a_radix_sorting_int_bather_merge/common/include/common.hpp"

namespace khruev_a_radix_sorting_int_bather_merge {

void KhruevARadixSortingIntBatherMergeSEQ::RadixSort(std::vector<int> &arr) {
  const int bits = 8;
  const int buckets = 1 << bits;
  const int mask = buckets - 1;
  const int passes = 32 / bits;

  std::vector<int> output(arr.size());

  for (int pass = 0; pass < passes; pass++) {
    std::vector<int> count(buckets, 0);

    int shift = pass * bits;

    for (int x : arr) {
      uint32_t ux = static_cast<uint32_t>(x) ^ 0x80000000U;
      int digit = (ux >> shift) & mask;
      count[digit]++;
    }

    for (int i = 1; i < buckets; i++) {
      count[i] += count[i - 1];
    }

    for (int i = arr.size() - 1; i >= 0; i--) {
      uint32_t ux = static_cast<uint32_t>(arr[i]) ^ 0x80000000U;
      int digit = (ux >> shift) & mask;
      output[--count[digit]] = arr[i];
    }

    arr = output;
  }
}

void KhruevARadixSortingIntBatherMergeSEQ::OddEvenMerge(std::vector<int> &a, int lo, int n, int r) {
  int step = r * 2;
  if (step < n) {
    OddEvenMerge(a, lo, n, step);
    OddEvenMerge(a, lo + r, n, step);

    for (int i = lo + r; i + r < lo + n; i += step) {
      if (a[i] > a[i + r]) {
        std::swap(a[i], a[i + r]);
      }
    }
  } else {
    if (a[lo] > a[lo + r]) {
      std::swap(a[lo], a[lo + r]);
    }
  }
}

void KhruevARadixSortingIntBatherMergeSEQ::OddEvenMergeSort(std::vector<int> &a, int lo, int n) {
  if (n > 1) {
    int m = n / 2;
    OddEvenMergeSort(a, lo, m);
    OddEvenMergeSort(a, lo + m, m);
    OddEvenMerge(a, lo, n, 1);
  }
}

KhruevARadixSortingIntBatherMergeSEQ::KhruevARadixSortingIntBatherMergeSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool KhruevARadixSortingIntBatherMergeSEQ::ValidationImpl() {
  return !GetInput().empty();
}

bool KhruevARadixSortingIntBatherMergeSEQ::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return true;
}

bool KhruevARadixSortingIntBatherMergeSEQ::RunImpl() {
  std::vector<int> data = GetInput();

  RadixSort(data);

  size_t original_size = data.size();

  size_t pow2 = 1;
  while (pow2 < original_size) {
    pow2 <<= 1;
  }

  data.resize(pow2, std::numeric_limits<int>::max());

  OddEvenMergeSort(data, 0, pow2);

  data.resize(original_size);

  GetOutput() = data;

  return true;
}

bool KhruevARadixSortingIntBatherMergeSEQ::PostProcessingImpl() {
  return GetOutput().size() == GetInput().size();
}

}  // namespace khruev_a_radix_sorting_int_bather_merge
