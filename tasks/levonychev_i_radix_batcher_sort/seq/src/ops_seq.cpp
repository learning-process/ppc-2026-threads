#include "levonychev_i_radix_batcher_sort/seq/include/ops_seq.hpp"

#include <cstddef>
#include <ranges>
#include <vector>

#include "levonychev_i_radix_batcher_sort/common/include/common.hpp"

namespace levonychev_i_radix_batcher_sort {

LevonychevIRadixBatcherSortSEQ::LevonychevIRadixBatcherSortSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

void LevonychevIRadixBatcherSortSEQ::CountingSort(InType &arr, size_t byte_index) {
  const size_t byte = 256;
  std::vector<int> count(byte, 0);
  OutType result(arr.size());

  bool is_last_byte = (byte_index == (sizeof(int) - 1ULL));
  for (auto number : arr) {
    int value_of_byte = (number >> (byte_index * 8ULL)) & 0xFF;

    if (is_last_byte) {
      value_of_byte ^= 0x80;
    }

    ++count[value_of_byte];
  }

  for (size_t i = 1ULL; i < byte; ++i) {
    count[i] += count[i - 1];
  }

  for (int &val : std::ranges::reverse_view(arr)) {
    int value_of_byte = (val >> (byte_index * 8ULL)) & 0xFF;

    if (is_last_byte) {
      value_of_byte ^= 0x80;
    }

    result[--count[value_of_byte]] = val;
  }
  arr = result;
}

bool LevonychevIRadixBatcherSortSEQ::ValidationImpl() {
  return !GetInput().empty();
}

bool LevonychevIRadixBatcherSortSEQ::PreProcessingImpl() {
  return true;
}

bool LevonychevIRadixBatcherSortSEQ::RunImpl() {
  GetOutput() = GetInput();

  for (size_t i = 0; i < sizeof(int); ++i) {
    CountingSort(GetOutput(), i);
  }

  return true;
}

bool LevonychevIRadixBatcherSortSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace levonychev_i_radix_batcher_sort
