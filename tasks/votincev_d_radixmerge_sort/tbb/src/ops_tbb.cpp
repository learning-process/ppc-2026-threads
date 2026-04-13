#include "votincev_d_radixmerge_sort/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <algorithm>
#include <array>
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

  uint32_t max_val = begin[0];
  for (int32_t i = 1; i < n; ++i) {
    max_val = std::max(begin[static_cast<size_t>(i)], max_val);
  }

  std::vector<uint32_t> buffer(static_cast<size_t>(n));
  uint32_t *src = begin;
  uint32_t *dst = buffer.data();

  //  int64_t для exp, чтобы избежать переполнения при exp * 10
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
void VotincevDRadixMergeSortTBB::Merge(uint32_t *data, int32_t left, int32_t mid, int32_t right, uint32_t *temp) {
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
  Merge(data, left, mid, right, temp);
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
