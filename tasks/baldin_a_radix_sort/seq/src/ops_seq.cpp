#include "baldin_a_radix_sort/seq/include/ops_seq.hpp"

#include <numeric>
#include <vector>

#include "baldin_a_radix_sort/common/include/common.hpp"
#include "util/include/util.hpp"

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
void countingSortByByte(std::vector<int> &arr, int byteIndex) {
  size_t n = arr.size();
  if (n == 0) {
    return;
  }

  std::vector<int> output(n);
  std::vector<int> count(256, 0);

  int shift = static_cast<int>(byteIndex) * 8;

  for (size_t i = 0; i < n; i++) {
    unsigned int rawVal = static_cast<unsigned int>(arr[i]);
    unsigned int byteVal = (rawVal >> shift) & 0xFF;

    if (byteIndex == sizeof(int) - 1) {
      byteVal ^= 128;
    }

    count[byteVal]++;
  }

  for (int i = 1; i < 256; i++) {
    count[i] += count[i - 1];
  }

  for (size_t i = n; i > 0; i--) {
    size_t idx = i - 1;
    unsigned int rawVal = static_cast<unsigned int>(arr[idx]);
    unsigned int byteVal = (rawVal >> shift) & 0xFF;

    if (byteIndex == sizeof(int) - 1) {
      byteVal ^= 128;
    }

    output[count[byteVal] - 1] = arr[idx];
    count[byteVal]--;
  }

  arr = output;
}
}  // namespace

bool BaldinARadixSortSEQ::RunImpl() {
  for (size_t byteIndex = 0; byteIndex < sizeof(int); byteIndex++) {
    countingSortByByte(GetOutput(), byteIndex);
  }
  return true;
}

bool BaldinARadixSortSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace baldin_a_radix_sort
