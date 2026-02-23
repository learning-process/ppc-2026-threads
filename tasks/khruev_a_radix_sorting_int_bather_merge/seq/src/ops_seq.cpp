#include "khruev_a_radix_sorting_int_bather_merge/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <vector>

#include "khruev_a_radix_sorting_int_bather_merge/common/include/common.hpp"

namespace khruev_a_radix_sorting_int_bather_merge {

void KhruevARadixSortingIntBatherMergeSEQ::compareExchange(std::vector<int> &a, int i, int j) {
  if (a[i] > a[j]) {
    std::swap(a[i], a[j]);
  }
}

void KhruevARadixSortingIntBatherMergeSEQ::oddEvenMerge(std::vector<int> &a, int lo, int n, int r) {
  int step = r * 2;
  if (step < n) {
    oddEvenMerge(a, lo, n, step);
    oddEvenMerge(a, lo + r, n, step);
    for (int i = lo + r; i + r < lo + n; i += step) {
      compareExchange(a, i, i + r);
    }
  } else {
    compareExchange(a, lo, lo + r);
  }
}

void KhruevARadixSortingIntBatherMergeSEQ::oddEvenMergeSort(std::vector<int> &a, int lo, int n) {
  if (n > 1) {
    int m = n / 2;
    oddEvenMergeSort(a, lo, m);
    oddEvenMergeSort(a, lo + m, m);
    oddEvenMerge(a, lo, n, 1);
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
  size_t n = GetInput().size();
  if (n == 0) {
    GetOutput().clear();
    return true;
  }

  size_t pow2 = 1;
  while (pow2 < n) {
    pow2 <<= 1;
  }

  std::vector<int> tmp = GetInput();
  tmp.resize(pow2, std::numeric_limits<int>::max());

  oddEvenMergeSort(tmp, 0, pow2);

  tmp.resize(n);
  GetOutput() = tmp;
  return true;
}

bool KhruevARadixSortingIntBatherMergeSEQ::PostProcessingImpl() {
  return GetOutput().size() == GetInput().size();
}

}  // namespace khruev_a_radix_sorting_int_bather_merge
