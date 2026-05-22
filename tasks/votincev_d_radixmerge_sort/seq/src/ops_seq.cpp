#include "votincev_d_radixmerge_sort/seq/include/ops_seq.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>
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
  std::array<int32_t, 10> count{};

  for (size_t i = 0; i < n; i++) {
    int32_t digit = (array[i] / exp) % 10;
    count.at(static_cast<size_t>(digit))++;
  }

  // префиксные суммы
  for (size_t i = 1; i < 10; i++) {
    count.at(i) += count.at(i - 1);
  }

  // формирую выходной массив
  for (int64_t i = static_cast<int64_t>(n) - 1; i >= 0; i--) {
    auto idx = static_cast<size_t>(i);
    auto digit = static_cast<size_t>((array.at(idx) / exp) % 10);
    size_t pos = static_cast<size_t>(count.at(digit)) - 1;
    output.at(pos) = array.at(idx);
    count.at(digit)--;
  }

  array = std::move(output);
}

bool VotincevDRadixMergeSortSEQ::RunImpl() {
  std::vector<int32_t> working_array = GetInput();

  auto [min_it, max_it] = std::ranges::minmax_element(working_array);
  int32_t min_val = *min_it;
  int32_t max_val = *max_it;

  // сдвиг в положительную область
  if (min_val < 0) {
    for (auto &num : working_array) {
      num -= min_val;
    }
    max_val -= min_val;
  }

  // цикл по разрядам
  for (int64_t exp = 1; static_cast<int64_t>(max_val) / exp > 0; exp *= 10) {
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
