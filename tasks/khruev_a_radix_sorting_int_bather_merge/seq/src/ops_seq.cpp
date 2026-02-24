#include "khruev_a_radix_sorting_int_bather_merge/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
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
      uint32_t digit = (ux >> shift) & mask;
      count[digit]++;
    }

    for (int i = 1; i < buckets; i++) {
      count[i] += count[i - 1];
    }

    for (size_t i = arr.size(); i-- > 0;) {
      uint32_t ux = static_cast<uint32_t>(arr[i]) ^ 0x80000000U;
      uint32_t digit = (ux >> shift) & mask;
      output[--count[digit]] = arr[i];
    }

    arr = output;
  }
}



void KhruevARadixSortingIntBatherMergeSEQ::OddEvenStage(
    std::vector<int>& a,
    size_t n,
    size_t p,
    size_t k) {

  const size_t block_size = 2 * p;

  for (size_t i = 0; i < n; ++i) {
    const size_t j = i ^ k;

    if (j > i && SameBlock(i, j, block_size)) {
      CompareExchange(a, i, j);
    }
  }
}

void KhruevARadixSortingIntBatherMergeSEQ::OddEvenMergeSort(std::vector<int>& a, size_t n) {
  for (size_t p = 1; p < n; p <<= 1) {
    for (size_t k = p; k > 0; k >>= 1) {
      OddEvenStage(a, n, p, k);
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

  RadixSort(data);

  size_t original_size = data.size();

  size_t pow2 = 1;
  while (pow2 < original_size) {
    pow2 <<= 1;
  }

  data.resize(pow2, std::numeric_limits<int>::max());

  OddEvenMergeSort(data, data.size());

  data.resize(original_size);

  GetOutput() = data;

  return true;
}

bool KhruevARadixSortingIntBatherMergeSEQ::PostProcessingImpl() {
  return GetOutput().size() == GetInput().size();
}

}  // namespace khruev_a_radix_sorting_int_bather_merge
