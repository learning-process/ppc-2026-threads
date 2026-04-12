#include "votincev_d_radixmerge_sort/stl/include/ops_stl.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <future>
#include <utility>
#include <vector>

#include "votincev_d_radixmerge_sort/common/include/common.hpp"

namespace votincev_d_radixmerge_sort {

VotincevDRadixMergeSortSTL::VotincevDRadixMergeSortSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool VotincevDRadixMergeSortSTL::ValidationImpl() {
  return !GetInput().empty();
}

bool VotincevDRadixMergeSortSTL::PreProcessingImpl() {
  return true;
}

// поразрядная сортировка для локальных блоков (LSD)
void VotincevDRadixMergeSortSTL::LocalRadixSort(uint32_t *begin, const uint32_t *end) {
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
      count[static_cast<size_t>((src[static_cast<size_t>(i)] / exp) % 10)]++;
    }
    for (int32_t i = 1; i < 10; ++i) {
      count[static_cast<size_t>(i)] += count[static_cast<size_t>(i - 1)];
    }
    for (int32_t i = n - 1; i >= 0; --i) {
      uint32_t digit = (src[static_cast<size_t>(i)] / exp) % 10;
      dst[--count[static_cast<size_t>(digit)]] = src[static_cast<size_t>(i)];
    }
    std::swap(src, dst);
  }

  if (src != begin) {
    std::copy(src, src + n, begin);
  }
}

// слияние двух отсортированных участков
void VotincevDRadixMergeSortSTL::Merge(uint32_t *data, int32_t left, int32_t mid, int32_t right, uint32_t *temp) {
  int32_t i = left;
  int32_t j = mid;
  int32_t k = left;
  while (i < mid && j < right) {
    temp[static_cast<size_t>(k++)] = (data[static_cast<size_t>(i)] <= data[static_cast<size_t>(j)])
                                         ? data[static_cast<size_t>(i++)]
                                         : data[static_cast<size_t>(j++)];
  }
  while (i < mid) {
    temp[static_cast<size_t>(k++)] = data[static_cast<size_t>(i++)];
  }
  while (j < right) {
    temp[static_cast<size_t>(k++)] = data[static_cast<size_t>(j++)];
  }

  std::copy(temp + left, temp + right, data + left);
}

// параллельная сортировка слиянием через std::async
void VotincevDRadixMergeSortSTL::ParallelRadixMergeSort(uint32_t *data, int32_t left, int32_t right,
                                                        uint32_t *temp) {  // NOLINT(misc-no-recursion)
  const int32_t grain_size = 4096;  // порог для перехода на последовательную сортировку

  if (right - left <= grain_size) {
    LocalRadixSort(data + left, data + right);
    return;
  }

  int32_t mid = left + ((right - left) / 2);

  // запускаем левую часть в отдельном потоке (аналог tbb::parallel_invoke)
  auto future = std::async(std::launch::async, [&] { ParallelRadixMergeSort(data, left, mid, temp); });

  // правую часть выполняем в текущем потоке
  ParallelRadixMergeSort(data, mid, right, temp);

  // ждем завершения левой части
  future.get();

  // сливаем результаты
  Merge(data, left, mid, right, temp);
}

bool VotincevDRadixMergeSortSTL::RunImpl() {
  const auto &input = GetInput();
  auto n = static_cast<int32_t>(input.size());

  // поиск минимума
  int32_t min_val = *std::ranges::min_element(input);

  // uint32_t чтобы избежать проблем с отрицательными числами
  std::vector<uint32_t> working_array(static_cast<size_t>(n));
  for (int32_t i = 0; i < n; ++i) {
    working_array[static_cast<size_t>(i)] =
        static_cast<uint32_t>(input[static_cast<size_t>(i)]) - static_cast<uint32_t>(min_val);
  }

  // параллельная сортировка
  std::vector<uint32_t> temp_buffer(static_cast<size_t>(n));
  ParallelRadixMergeSort(working_array.data(), 0, n, temp_buffer.data());

  // восстановление исходных значений
  std::vector<int32_t> result(static_cast<size_t>(n));
  for (int32_t i = 0; i < n; ++i) {
    result[static_cast<size_t>(i)] =
        static_cast<int32_t>(working_array[static_cast<size_t>(i)] + static_cast<uint32_t>(min_val));
  }

  GetOutput() = std::move(result);
  return true;
}

bool VotincevDRadixMergeSortSTL::PostProcessingImpl() {
  return true;
}

}  // namespace votincev_d_radixmerge_sort
