#include "baldin_a_radix_sort/seq/include/ops_seq.hpp"

#include <cstddef>
#include <vector>

#include "baldin_a_radix_sort/common/include/common.hpp"

namespace baldin_a_radix_sort {

BaldinARadixSortSEQ::BaldinARadixSortSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool BaldinARadixSortSEQ::ValidationImpl() {
  return true;
}

bool BaldinARadixSortSEQ::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

namespace {
void CountingSortByByte(std::vector<int> &arr, size_t byte_index) {
  size_t n = arr.size();
  if (n == 0) {
    return;
  }

  std::vector<int> output(n);
  std::vector<int> count(256, 0);

  size_t shift = byte_index * 8;

  for (size_t i = 0; i < n; i++) {
    auto raw_val = static_cast<unsigned int>(arr[i]);
    unsigned int byte_val = (raw_val >> shift) & 0xFF;

    if (byte_index == sizeof(int) - 1) {
      byte_val ^= 128;
    }

    count[byte_val]++;
  }

  for (int i = 1; i < 256; i++) {
    count[i] += count[i - 1];
  }

  for (size_t i = n; i > 0; i--) {
    size_t idx = i - 1;
    auto raw_val = static_cast<unsigned int>(arr[idx]);
    unsigned int byte_val = (raw_val >> shift) & 0xFF;

    if (byte_index == sizeof(int) - 1) {
      byte_val ^= 128;
    }

    output[count[byte_val] - 1] = arr[idx];
    count[byte_val]--;
  }

  arr = output;
}
}  // namespace

bool BaldinARadixSortSEQ::RunImpl() {
  for (size_t byte_index = 0; byte_index < sizeof(int); byte_index++) {
    CountingSortByByte(GetOutput(), byte_index);
  }
  return true;
}

bool BaldinARadixSortSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace baldin_a_radix_sort
