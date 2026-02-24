#include "khruev_a_radix_sorting_int_bather_merge/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <vector>

#include "khruev_a_radix_sorting_int_bather_merge/common/include/common.hpp"

namespace khruev_a_radix_sorting_int_bather_merge {

void KhruevARadixSortingIntBatherMergeSEQ::radixSort(std::vector<int> &arr) {
  const int BITS = 8;
  const int BUCKETS = 1 << BITS;  // 256
  const int MASK = BUCKETS - 1;
  const int PASSES = 32 / BITS;  // 4 прохода

  std::vector<int> output(arr.size());

  for (int pass = 0; pass < PASSES; pass++) {
    std::vector<int> count(BUCKETS, 0);

    int shift = pass * BITS;

    // Подсчёт
    for (int x : arr) {
      uint32_t ux = static_cast<uint32_t>(x) ^ 0x80000000u;
      int digit = (ux >> shift) & MASK;
      count[digit]++;
    }

    // Префиксные суммы
    for (int i = 1; i < BUCKETS; i++) {
      count[i] += count[i - 1];
    }

    // Стабильная сортировка
    for (int i = arr.size() - 1; i >= 0; i--) {
      uint32_t ux = static_cast<uint32_t>(arr[i]) ^ 0x80000000u;
      int digit = (ux >> shift) & MASK;
      output[--count[digit]] = arr[i];
    }

    arr = output;
  }
}

void KhruevARadixSortingIntBatherMergeSEQ::oddEvenMerge(std::vector<int> &a, int lo, int n, int r) {
  int step = r * 2;
  if (step < n) {
    oddEvenMerge(a, lo, n, step);
    oddEvenMerge(a, lo + r, n, step);

    for (int i = lo + r; i + r < lo + n; i += step) {
      if (a[i] > a[i + r]) {
        std::swap(a[i], a[i + r]);
      }
    }
  } else {
    if (a[lo] > a[lo + r]) {
      std::swap(a[lo], a[lo + r]);
    }
  }
}

void KhruevARadixSortingIntBatherMergeSEQ::oddEvenMergeSort(std::vector<int> &a, int lo, int n) {
  if (n > 1) {
    int m = n / 2;
    oddEvenMergeSort(a, lo, m);
    oddEvenMergeSort(a, lo + m, m);
    oddEvenMerge(a, lo, n, 1);
  }
}

KhruevARadixSortingIntBatherMergeSEQ::KhruevARadixSortingIntBatherMergeSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool KhruevARadixSortingIntBatherMergeSEQ::ValidationImpl() {
  return !GetInput().empty();
}

bool KhruevARadixSortingIntBatherMergeSEQ::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return true;
}

bool KhruevARadixSortingIntBatherMergeSEQ::RunImpl() {
  // 1️⃣ Копируем вход в локальный вектор
  std::vector<int> data = GetInput();

  // 2️⃣ Сначала локальная сортировка radix
  radixSort(data);

  // 3️⃣ Запоминаем реальный размер
  size_t original_size = data.size();

  // 4️⃣ Делаем размер степенью двойки
  size_t pow2 = 1;
  while (pow2 < original_size) {
    pow2 <<= 1;
  }

  // 5️⃣ Дополняем "бесконечностями"
  data.resize(pow2, std::numeric_limits<int>::max());

  // 6️⃣ Теперь можно безопасно вызывать Batcher
  oddEvenMergeSort(data, 0, pow2);

  // 7️⃣ Убираем добавленные элементы
  data.resize(original_size);

  // 8️⃣ Кладём результат в output
  GetOutput() = data;

  return true;
}

bool KhruevARadixSortingIntBatherMergeSEQ::PostProcessingImpl() {
  return GetOutput().size() == GetInput().size();
}

}  // namespace khruev_a_radix_sorting_int_bather_merge
