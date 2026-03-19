#include "khruev_a_radix_sorting_int_bather_merge/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include "khruev_a_radix_sorting_int_bather_merge/common/include/common.hpp"

namespace khruev_a_radix_sorting_int_bather_merge {

void KhruevARadixSortingIntBatherMergeSEQ::CompareExchange(std::vector<int> &a, size_t i, size_t j) {
  if (a[i] > a[j]) {
    std::swap(a[i], a[j]);
  }
}

void KhruevARadixSortingIntBatherMergeSEQ::RadixSort(std::vector<int> &arr) {
  const int bits = 8;
  const int buckets = 1 << bits;
  const int mask = buckets - 1;
  const int passes = 32 / bits;

  std::vector<int> temp(arr.size());
  std::vector<int> *src = &arr;
  std::vector<int> *dst = &temp;

  for (int pass = 0; pass < passes; pass++) {
    std::vector<int> count(buckets, 0);
    int shift = pass * bits;

    for (int x : *src) {
      uint32_t ux = static_cast<uint32_t>(x) ^ 0x80000000U;
      uint32_t digit = (ux >> shift) & mask;
      count[digit]++;
    }

    for (int i = 1; i < buckets; i++) {
      count[i] += count[i - 1];
    }

    for (size_t i = src->size(); i-- > 0;) {
      uint32_t ux = static_cast<uint32_t>((*src)[i]) ^ 0x80000000U;
      uint32_t digit = (ux >> shift) & mask;
      (*dst)[--count[digit]] = (*src)[i];
    }

    std::swap(src, dst);
  }
}

void KhruevARadixSortingIntBatherMergeSEQ::OddEvenMerge(std::vector<int> &a, size_t n) {
  for (size_t po = n / 2; po > 0; po >>= 1) {
    if (po == n / 2) {
      for (size_t i = 0; i < po; ++i) {
        CompareExchange(a, i, i + po);
      }
    } else {
      for (size_t i = po; i < n - po; i += 2 * po) {
        for (size_t j = 0; j < po; ++j) {
          CompareExchange(a, i + j, i + j + po);
        }
      }
    }
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
  size_t original_size = data.size();

  if (original_size <= 1) {
    GetOutput() = data;
    return true;
  }

  size_t pow2 = 1;
  while (pow2 < original_size) {
    pow2 <<= 1;
  }

  data.resize(pow2, std::numeric_limits<int>::max());

  size_t half = pow2 / 2;

  auto half_dist = static_cast<std::ptrdiff_t>(half);
  std::vector<int> left(data.begin(), data.begin() + half_dist);
  std::vector<int> right(data.begin() + half_dist, data.end());

  RadixSort(left);
  RadixSort(right);

  std::ranges::copy(left, data.begin());
  std::ranges::copy(right, data.begin() + half_dist);

  OddEvenMerge(data, data.size());

  data.resize(original_size);

  GetOutput() = data;

  return true;
}

bool KhruevARadixSortingIntBatherMergeSEQ::PostProcessingImpl() {
  return GetOutput().size() == GetInput().size();
}

}  // namespace khruev_a_radix_sorting_int_bather_merge
