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

// проверка входных данных
bool VotincevDRadixMergeSortSEQ::ValidationImpl() {
  // проверка: входной вектор не должен быть пустым
  return !GetInput().empty();
}

// препроцессинг
bool VotincevDRadixMergeSortSEQ::PreProcessingImpl() {
  return true;
}

// вспомогательный метод для распределения и слияния разрядов
void VotincevDRadixMergeSortSEQ::SortByDigit(std::vector<int32_t> &array, int32_t exp) {
  std::vector<std::vector<int32_t>> buckets(10);

  // распределение элементов по корзинам
  for (const auto &num : array) {
    int32_t digit = (num / exp) % 10;
    buckets[digit].push_back(num);
  }

  // простое слияние корзин обратно в рабочий массив
  size_t index = 0;
  for (int i = 0; i < 10; ++i) {
    for (const auto &val : buckets[i]) {
      array[index++] = val;
    }
    // очистка корзины для следующего разряда
    buckets[i].clear();
  }
}

// основной метод алгоритма
bool VotincevDRadixMergeSortSEQ::RunImpl() {
  if (GetInput().empty()) {
    return false;
  }

  // локальная копия данных для сортировки
  std::vector<int32_t> working_array = GetInput();

  // обработка отрицательных чисел
  int32_t min_val = *std::ranges::min_element(working_array);

  if (min_val < 0) {
    for (auto &num : working_array) {
      num -= min_val;
    }
  }

  // ищем максимальное число для определения количества разрядов
  int32_t max_val = *std::ranges::max_element(working_array);

  // цикл по разрядам (единицы, десятки, сотни...)
  for (int32_t exp = 1; max_val / exp > 0; exp *= 10) {
    SortByDigit(working_array, exp);
  }

  // возвращаем значения к исходному диапазону
  if (min_val < 0) {
    for (auto &num : working_array) {
      num += min_val;
    }
  }

  // запись результата в выходные данные
  GetOutput() = working_array;

  return true;
}

// постпроцессинг
bool VotincevDRadixMergeSortSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace votincev_d_radixmerge_sort
