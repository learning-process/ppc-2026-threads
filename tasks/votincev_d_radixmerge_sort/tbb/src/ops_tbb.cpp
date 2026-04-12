#include "votincev_d_radixmerge_sort/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <algorithm>
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
  int32_t n = static_cast<int32_t>(end - begin);
  if (n <= 1) {
    return;
  }

  uint32_t max_val = begin[0];
  for (int32_t i = 1; i < n; ++i) {
    if (begin[i] > max_val) {
      max_val = begin[i];
    }
  }

  std::vector<uint32_t> buffer(n);
  uint32_t *src = begin;
  uint32_t *dst = buffer.data();

  //  int64_t для exp, чтобы избежать переполнения при exp * 10
  for (int64_t exp = 1; static_cast<int64_t>(max_val) / exp > 0; exp *= 10) {
    int32_t count[10] = {0};

    for (int32_t i = 0; i < n; ++i) {
      count[(src[i] / exp) % 10]++;
    }
    for (int32_t i = 1; i < 10; ++i) {
      count[i] += count[i - 1];
    }
    for (int32_t i = n - 1; i >= 0; --i) {
      uint32_t digit = (src[i] / exp) % 10;
      dst[--count[digit]] = src[i];
    }
    std::swap(src, dst);
  }

  if (src != begin) {
    std::copy(src, src + n, begin);
  }
}

// слияние двух отсортированных участков
void VotincevDRadixMergeSortTBB::Merge(uint32_t *data, int32_t left, int32_t mid, int32_t right, uint32_t *temp) {
  int32_t i = left, j = mid, k = left;
  while (i < mid && j < right) {
    temp[k++] = (data[i] <= data[j]) ? data[i++] : data[j++];
  }
  while (i < mid) {
    temp[k++] = data[i++];
  }
  while (j < right) {
    temp[k++] = data[j++];
  }

  std::copy(temp + left, temp + right, data + left);
}

// параллельная сортировка слиянием(сортировка + слияние)
void VotincevDRadixMergeSortTBB::ParallelRadixMergeSort(uint32_t *data, int32_t left, int32_t right, uint32_t *temp) {
  const int32_t GRAIN_SIZE = 4096;  // порог для перехода на последовательную сортировку

  if (right - left <= GRAIN_SIZE) {
    LocalRadixSort(data + left, data + right);
    return;
  }

  int32_t mid = left + (right - left) / 2;

  // рекурсивно запускаются две задачи в параллель
  tbb::parallel_invoke([&] { ParallelRadixMergeSort(data, left, mid, temp); },
                       [&] { ParallelRadixMergeSort(data, mid, right, temp); });

  // слияние результатов
  Merge(data, left, mid, right, temp);
}

bool VotincevDRadixMergeSortTBB::RunImpl() {
  const auto &input = GetInput();
  int32_t n = static_cast<int32_t>(input.size());

  // поиск минимума
  int32_t min_val = tbb::parallel_reduce(tbb::blocked_range<int32_t>(0, n), input[0],
                                         [&](const tbb::blocked_range<int32_t> &r, int32_t local_min) {
    for (int32_t i = r.begin(); i < r.end(); ++i) {
      if (input[i] < local_min) {
        local_min = input[i];
      }
    }
    return local_min;
  }, [](int32_t a, int32_t b) { return std::min(a, b); });

  // приведение к положительным uint32
  std::vector<uint32_t> working_array(n);
  tbb::parallel_for(
      0, n, [&](int32_t i) { working_array[i] = static_cast<uint32_t>(input[i]) - static_cast<uint32_t>(min_val); });

  // параллельная сортировка со слиянием
  std::vector<uint32_t> temp_buffer(n);
  ParallelRadixMergeSort(working_array.data(), 0, n, temp_buffer.data());

  // восстановление исходных значений
  std::vector<int32_t> result(n);
  tbb::parallel_for(
      0, n, [&](int32_t i) { result[i] = static_cast<int32_t>(working_array[i] + static_cast<uint32_t>(min_val)); });

  GetOutput() = std::move(result);
  return true;
}

bool VotincevDRadixMergeSortTBB::PostProcessingImpl() {
  return true;
}

}  // namespace votincev_d_radixmerge_sort
