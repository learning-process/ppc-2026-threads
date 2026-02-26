#include "frolova_s_radix_sort_double/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstring>
#include <vector>

namespace frolova_s_radix_sort_double {

FrolovaSRadixSortDoubleSEQ::FrolovaSRadixSortDoubleSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;  // сохраняем входной вектор
}

bool FrolovaSRadixSortDoubleSEQ::ValidationImpl() {
  // Проверяем, что входной вектор не пуст
  return !GetInput().empty();
}

bool FrolovaSRadixSortDoubleSEQ::PreProcessingImpl() {
  // Дополнительная подготовка не требуется
  return true;
}

bool FrolovaSRadixSortDoubleSEQ::RunImpl() {
  const std::vector<double>& input = GetInput();
  if (input.empty()) return false;

  std::vector<double> working = input;

  // Параметры поразрядной сортировки (по байтам)
  const int radix = 256;          // значения байта 0..255
  const int num_bits = 8;
  const int num_passes = sizeof(uint64_t); // 8 проходов для double

  std::vector<int> count(radix);
  std::vector<double> temp(working.size());

  for (int pass = 0; pass < num_passes; ++pass) {
    std::fill(count.begin(), count.end(), 0);

    // Подсчёт количества элементов для каждого значения текущего байта
    for (double val : working) {
      uint64_t bits;
      std::memcpy(&bits, &val, sizeof(double));
      int byte = (bits >> (pass * num_bits)) & 0xFF;
      ++count[byte];
    }

    // Преобразование счётчиков в начальные позиции
    int total = 0;
    for (int i = 0; i < radix; ++i) {
      int old = count[i];
      count[i] = total;
      total += old;
    }

    // Распределение элементов во временный массив
    for (double val : working) {
      uint64_t bits;
      std::memcpy(&bits, &val, sizeof(double));
      int byte = (bits >> (pass * num_bits)) & 0xFF;
      temp[count[byte]++] = val;
    }

    // Копирование обратно
    working = temp;
  }

  // Коррекция порядка отрицательных чисел
  std::vector<double> negative, positive;
  for (double val : working) {
    if (val < 0)
      negative.push_back(val);
    else
      positive.push_back(val);
  }
  std::reverse(negative.begin(), negative.end());

  // Простое слияние: сначала отрицательные, потом положительные
  working.clear();
  working.insert(working.end(), negative.begin(), negative.end());
  working.insert(working.end(), positive.begin(), positive.end());

  // Сохраняем результат
  GetOutput() = std::move(working);
  return true;
}

bool FrolovaSRadixSortDoubleSEQ::PostProcessingImpl() {
  // Результат уже сохранён в GetOutput()
  return true;
}

}  // namespace frolova_s_radix_sort_double