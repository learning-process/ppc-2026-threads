#include "levonychev_i_radix_batcher_sort/seq/include/ops_seq.hpp"

#include <numeric>
#include <vector>

#include "levonychev_i_radix_batcher_sort/common/include/common.hpp"
#include "util/include/util.hpp"

namespace levonychev_i_radix_batcher_sort {

LevonychevIRadixBatcherSortSEQ::LevonychevIRadixBatcherSortSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

void LevonychevIRadixBatcherSortSEQ::CountingSort(InType &arr, int byte_index) {
  const int byte = 256;
  std::vector<int> count(byte, 0);
  OutType result(arr.size());

  bool is_last_byte = (byte_index == (sizeof(int) - 1));
  for (auto number : arr) {
    int value_of_byte = (number >> (byte_index * 8)) & 0xFF;

    if (is_last_byte) {
      value_of_byte ^= 0x80;
    }

    ++count[value_of_byte];
  }

  for (int i = 1; i < byte; ++i) {
    count[i] += count[i - 1];
  }

  for (int i = arr.size() - 1; i >= 0; --i) {
    int value_of_byte = (arr[i] >> (byte_index * 8)) & 0xFF;

    if (is_last_byte) {
      value_of_byte ^= 0x80;
    }

    result[--count[value_of_byte]] = arr[i];
  }
  arr = result;
}

bool LevonychevIRadixBatcherSortSEQ::ValidationImpl() {
  return true;
}

bool LevonychevIRadixBatcherSortSEQ::PreProcessingImpl() {
  return true;
}

bool LevonychevIRadixBatcherSortSEQ::RunImpl() {
  GetOutput() = GetInput();

  for (const auto &i : GetOutput()) {
    std::cout << i << ' ';
  }
  std::cout << std::endl;

  for (int32_t i = 0; i < sizeof(int); ++i) {
    CountingSort(GetOutput(), i);
  }

  for (const auto &i : GetOutput()) {
    std::cout << i << ' ';
  }
  std::cout << std::endl;

  return true;
}

bool LevonychevIRadixBatcherSortSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace levonychev_i_radix_batcher_sort
