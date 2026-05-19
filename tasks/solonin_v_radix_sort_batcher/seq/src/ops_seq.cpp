#include "solonin_v_radix_sort_batcher/seq/include/ops_seq.hpp"
#include "solonin_v_radix_sort_batcher/common/include/common.hpp"
#include <cstddef>
#include <vector>

namespace solonin_v_radix_sort_batcher {

RadixSortBatcherSEQ::RadixSortBatcherSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

void RadixSortBatcherSEQ::SortByDigit(std::vector<int> &data, size_t pos) {
  const size_t k_base = 256;
  std::vector<int> freq(k_base, 0);
  std::vector<int> out(data.size());
  bool last = (pos == sizeof(int) - 1ULL);

  for (int v : data) {
    int byte_val = (v >> (pos * 8ULL)) & 0xFF;
    if (last) {
      byte_val ^= 0x80;
    }
    freq[byte_val]++;
  }
  for (size_t i = 1; i < k_base; ++i) {
    freq[i] += freq[i - 1];
  }
  for (int i = static_cast<int>(data.size()) - 1; i >= 0; --i) {
    int byte_val = (data[i] >> (pos * 8ULL)) & 0xFF;
    if (last) {
      byte_val ^= 0x80;
    }
    out[--freq[byte_val]] = data[i];
  }
  data = out;
}

void RadixSortBatcherSEQ::RadixSort(std::vector<int> &data) {
  for (size_t pos_idx = 0; pos_idx < sizeof(int); ++pos_idx) {
    SortByDigit(data, pos_idx);
  }
}

bool RadixSortBatcherSEQ::ValidationImpl() {
  return !GetInput().empty();
}

bool RadixSortBatcherSEQ::PreProcessingImpl() {
  return true;
}

bool RadixSortBatcherSEQ::RunImpl() {
  GetOutput() = GetInput();
  RadixSort(GetOutput());
  return true;
}

bool RadixSortBatcherSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace solonin_v_radix_sort_batcher
