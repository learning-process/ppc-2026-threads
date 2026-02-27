#include "spichek_d_radix_sort_for_integers_with_simple_merging/seq/include/ops_seq.hpp"

#include <algorithm>
#include <vector>

namespace spichek_d_radix_sort_for_integers_with_simple_merging {

RadixSortSEQ::RadixSortSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool RadixSortSEQ::ValidationImpl() {
  // Проверяем, что входные данные не пусты (или любые другие условия)
  return true;
}

bool RadixSortSEQ::PreProcessingImpl() {
  // Копируем входные данные в выходной буфер для последующей сортировки
  GetOutput() = GetInput();
  return true;
}

bool RadixSortSEQ::RunImpl() {
  if (GetOutput().empty()) {
    return true;
  }

  RadixSort(GetOutput());
  return true;
}

bool RadixSortSEQ::PostProcessingImpl() {
  // Проверка: отсортирован ли массив
  return std::is_sorted(GetOutput().begin(), GetOutput().end());
}

void RadixSortSEQ::RadixSort(std::vector<int> &data) {
  if (data.empty()) {
    return;
  }

  // Обработка отрицательных чисел: найдем минимум и сместим все значения,
  // чтобы работать только с положительными числами (LSD Radix Sort)
  int min_val = *std::min_element(data.begin(), data.end());
  if (min_val < 0) {
    for (auto &x : data) {
      x -= min_val;
    }
  }

  int max_val = *std::max_element(data.begin(), data.end());

  // Сама поразрядная сортировка (LSD) по основанию 10
  for (int exp = 1; max_val / exp > 0; exp *= 10) {
    std::vector<int> output(data.size());
    int count[10] = {0};

    for (int x : data) {
      count[(x / exp) % 10]++;
    }

    for (int i = 1; i < 10; i++) {
      count[i] += count[i - 1];
    }

    for (int i = static_cast<int>(data.size()) - 1; i >= 0; i--) {
      output[count[(data[i] / exp) % 10] - 1] = data[i];
      count[(data[i] / exp) % 10]--;
    }

    data = output;
  }

  // Возвращаем смещение назад, если оно было
  if (min_val < 0) {
    for (auto &x : data) {
      x += min_val;
    }
  }
}

}  // namespace spichek_d_radix_sort_for_integers_with_simple_merging
