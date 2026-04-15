#include "votincev_d_radixmerge_sort/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "votincev_d_radixmerge_sort/common/include/common.hpp"

namespace votincev_d_radixmerge_sort {

VotincevDRadixMergeSortTBB::VotincevDRadixMergeSortTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool VotincevDRadixMergeSortTBB::ValidationImpl() {
  return !GetInput().empty();
}

bool VotincevDRadixMergeSortTBB::PreProcessingImpl() {
  return true;
}

// поразрядная сортировка для локальных блоков
void VotincevDRadixMergeSortTBB::LocalRadixSort(uint32_t *begin, uint32_t *end) {
  auto n = static_cast<int32_t>(end - begin);
  if (n <= 1) {
    return;
  }

  uint32_t max_val = *std::max_element(begin, end);

  std::vector<uint32_t> buffer(static_cast<size_t>(n));
  uint32_t *src = begin;
  uint32_t *dst = buffer.data();

  for (int64_t exp = 1; (static_cast<int64_t>(max_val) / exp) > 0; exp *= 10) {
    std::array<int32_t, 10> count{};

    for (int32_t i = 0; i < n; ++i) {
      count.at(static_cast<size_t>((src[i] / exp) % 10))++;
    }
    for (int32_t i = 1; i < 10; ++i) {
      count.at(static_cast<size_t>(i)) += count.at(static_cast<size_t>(i - 1));
    }
    for (int32_t i = n - 1; i >= 0; --i) {
      auto digit = static_cast<size_t>((src[i] / exp) % 10);

      size_t target_idx = static_cast<size_t>(count.at(digit)) - 1;
      dst[target_idx] = src[i];

      count.at(digit)--;
    }
    std::swap(src, dst);
  }

  if (src != begin) {
    std::copy(src, src + n, begin);
  }
}

// слияние двух отсортированных участков
void VotincevDRadixMergeSortTBB::Merge(const uint32_t *src, uint32_t *dst, int32_t left, int32_t mid, int32_t right) {
  int32_t i = left;
  int32_t j = mid;
  int32_t k = left;
  while (i < mid && j < right) {
    dst[k++] = (src[i] <= src[j]) ? src[i++] : src[j++];
  }
  while (i < mid) {
    dst[k++] = src[i++];
  }
  while (j < right) {
    dst[k++] = src[j++];
  }
}

// параллельная сортировка слиянием(сортировка + слияние)
void VotincevDRadixMergeSortTBB::ParallelRadixMergeSort(uint32_t *data, int32_t left, int32_t right, uint32_t *temp) {
  const int32_t grain_size = 4096;  // порог для перехода на последовательную сортировку

  if (right - left <= grain_size) {
    LocalRadixSort(data + left, data + right);
    return;
  }

  int32_t mid = left + ((right - left) / 2);

  // рекурсивно запускаются две задачи в параллель
  tbb::parallel_invoke([&] { ParallelRadixMergeSort(data, left, mid, temp); },
                       [&] { ParallelRadixMergeSort(data, mid, right, temp); });

  // слияние результатов
  Merge(data, temp, left, mid, right);

  // копируем в основной массив
  std::copy(temp + left, temp + right, data + left);
}

bool VotincevDRadixMergeSortTBB::RunImpl() {
  const auto &input = GetInput();
  auto n = static_cast<int32_t>(input.size());

  // поиск минимума
  int32_t min_val = tbb::parallel_reduce(tbb::blocked_range<int32_t>(0, n), input[0],
                                         [&](const tbb::blocked_range<int32_t> &r, int32_t local_min) {
    for (int32_t i = r.begin(); i < r.end(); ++i) {
      local_min = std::min(input[static_cast<size_t>(i)], local_min);
    }
    return local_min;
  }, [](int32_t a, int32_t b) { return std::min(a, b); });

  // приведение к положительным uint32
  std::vector<uint32_t> working_array(static_cast<size_t>(n));
  tbb::parallel_for(0, n, [&](int32_t i) {
    working_array[static_cast<size_t>(i)] =
        static_cast<uint32_t>(input[static_cast<size_t>(i)]) - static_cast<uint32_t>(min_val);
  });

  // параллельная сортировка со слиянием
  std::vector<uint32_t> temp_buffer(static_cast<size_t>(n));
  ParallelRadixMergeSort(working_array.data(), 0, n, temp_buffer.data());

  // восстановление исходных значений
  std::vector<int32_t> result(static_cast<size_t>(n));
  tbb::parallel_for(0, n, [&](int32_t i) {
    result[static_cast<size_t>(i)] =
        static_cast<int32_t>(working_array[static_cast<size_t>(i)] + static_cast<uint32_t>(min_val));
  });

  GetOutput() = std::move(result);
  return true;
}

bool VotincevDRadixMergeSortTBB::PostProcessingImpl() {
  return true;
}

}  // namespace votincev_d_radixmerge_sort
