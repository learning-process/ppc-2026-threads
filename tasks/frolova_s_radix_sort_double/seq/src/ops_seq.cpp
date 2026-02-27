#include "frolova_s_radix_sort_double/seq/include/ops_seq.hpp"

#include <algorithm>
#include <bit>      // для std::bit_cast (C++20)
#include <cstdint>  // для uint64_t
#include <cstring>  // может остаться, если не используем bit_cast
#include <utility>  // для std::move
#include <vector>

namespace frolova_s_radix_sort_double {

FrolovaSRadixSortDoubleSEQ::FrolovaSRadixSortDoubleSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool FrolovaSRadixSortDoubleSEQ::ValidationImpl() {
  return !GetInput().empty();
}

bool FrolovaSRadixSortDoubleSEQ::PreProcessingImpl() {
  return true;
}

bool FrolovaSRadixSortDoubleSEQ::RunImpl() {
  const std::vector<double> &input = GetInput();
  if (input.empty()) {
    return false;
  }

  std::vector<double> working = input;

  const int radix = 256;
  const int num_bits = 8;
  const int num_passes = sizeof(uint64_t);  // 8 проходов

  std::vector<int> count(radix);
  std::vector<double> temp(working.size());

  for (int pass = 0; pass < num_passes; ++pass) {
    std::ranges::fill(count, 0);  // вместо std::fill

    // Подсчёт
    for (double val : working) {
      uint64_t bits = std::bit_cast<uint64_t>(val);  // безопасное преобразование
      int byte = static_cast<int>((bits >> (pass * num_bits)) & 0xFF);
      ++count[byte];
    }

    // Преобразование счётчиков в позиции
    int total = 0;
    for (int i = 0; i < radix; ++i) {
      int old = count[i];
      count[i] = total;
      total += old;
    }

    // Распределение
    for (double val : working) {
      uint64_t bits = std::bit_cast<uint64_t>(val);
      int byte = static_cast<int>((bits >> (pass * num_bits)) & 0xFF);
      temp[count[byte]++] = val;
    }

    working = temp;
  }

  // Коррекция порядка отрицательных чисел
  std::vector<double> negative;
  std::vector<double> positive;
  for (double val : working) {
    if (val < 0) {
      negative.push_back(val);
    } else {
      positive.push_back(val);
    }
  }
  std::ranges::reverse(negative);  // вместо std::reverse

  working.clear();
  working.insert(working.end(), negative.begin(), negative.end());
  working.insert(working.end(), positive.begin(), positive.end());

  GetOutput() = std::move(working);
  return true;
}

bool FrolovaSRadixSortDoubleSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace frolova_s_radix_sort_double
