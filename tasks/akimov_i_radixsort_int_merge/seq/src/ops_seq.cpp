#include "akimov_i_radixsort_int_merge/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstdint>
#include <vector>

#include "akimov_i_radixsort_int_merge/common/include/common.hpp"
#include "util/include/util.hpp"

namespace akimov_i_radixsort_int_merge {

AkimovIRadixSortIntMergeSEQ::AkimovIRadixSortIntMergeSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType();
}

bool AkimovIRadixSortIntMergeSEQ::ValidationImpl() {
  return !GetInput().empty();
}

bool AkimovIRadixSortIntMergeSEQ::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool AkimovIRadixSortIntMergeSEQ::RunImpl() {
  auto &arr = GetOutput();
  if (arr.empty()) return true;

  constexpr uint32_t sign_mask = 0x80000000;
  for (int &x : arr) {
    x ^= sign_mask;
  }

  const int num_bytes = static_cast<int>(sizeof(int));
  std::vector<int> temp(arr.size());

  for (int byte_pos = 0; byte_pos < num_bytes; ++byte_pos) {
    int count[256] = {0};
    for (int x : arr) {
      uint8_t byte = (x >> (byte_pos * 8)) & 0xFF;
      ++count[byte];
    }

    for (int i = 1; i < 256; ++i) {
      count[i] += count[i - 1];
    }

    for (int i = static_cast<int>(arr.size()) - 1; i >= 0; --i) {
      uint8_t byte = (arr[i] >> (byte_pos * 8)) & 0xFF;
      temp[--count[byte]] = arr[i];
    }

    arr.swap(temp);
  }

  for (int &x : arr) {
    x ^= sign_mask;
  }

  return true;
}

bool AkimovIRadixSortIntMergeSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace akimov_i_radixsort_int_merge
