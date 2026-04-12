#include "votincev_d_radixmerge_sort/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "votincev_d_radixmerge_sort/common/include/common.hpp"

namespace votincev_d_radixmerge_sort {

VotincevDRadixMergeSortSEQ::VotincevDRadixMergeSortSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool VotincevDRadixMergeSortSEQ::ValidationImpl() {
  return !GetInput().empty();
}

bool VotincevDRadixMergeSortSEQ::PreProcessingImpl() {
  return true;
}

// поразрядная сортировка
void VotincevDRadixMergeSortSEQ::SortByDigit(std::vector<int32_t> &array, int32_t exp) {
  size_t n = array.size();
  std::vector<int32_t> output(n);
  int32_t count[10] = {0};

  for (size_t i = 0; i < n; i++) {
    int32_t digit = (array[i] / exp) % 10;
    count[digit]++;
  }

  // префиксные суммы
  for (int i = 1; i < 10; i++) {
    count[i] += count[i - 1];
  }

  // теперь count[i] содержит позицию, перед которой заканчиваются элементы с цифрой i

  // формирую выходной массив
  for (int64_t i = n - 1; i >= 0; i--) {
    int32_t digit = (array[i] / exp) % 10;
    output[count[digit] - 1] = array[i];
    count[digit]--;
  }

  array = std::move(output);
}

bool VotincevDRadixMergeSortSEQ::RunImpl() {
  std::vector<int32_t> working_array = GetInput();

  // поиск min и max за один проход
  auto [min_it, max_it] = std::minmax_element(working_array.begin(), working_array.end());
  int32_t min_val = *min_it;
  int32_t max_val = *max_it;

  // сдвиг в положительную область
  if (min_val < 0) {
    for (auto &num : working_array) {
      num -= min_val;
    }
    max_val -= min_val;
  }

  // цикл по разрядам, int64_t для exp, чтобы избежать переполнения при exp * 10
  for (int64_t exp = 1; max_val / exp > 0; exp *= 10) {
    SortByDigit(working_array, static_cast<int32_t>(exp));
  }

  // возврат к исходному диапазону
  if (min_val < 0) {
    for (auto &num : working_array) {
      num += min_val;
    }
  }

  GetOutput() = std::move(working_array);
  return true;
}

bool VotincevDRadixMergeSortSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace votincev_d_radixmerge_sort
