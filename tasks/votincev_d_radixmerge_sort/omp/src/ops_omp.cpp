#include "votincev_d_radixmerge_sort/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cstring>
#include <vector>

namespace votincev_d_radixmerge_sort {

VotincevDRadixMergeSortOMP::VotincevDRadixMergeSortOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool VotincevDRadixMergeSortOMP::ValidationImpl() {
  return !GetInput().empty();
}

bool VotincevDRadixMergeSortOMP::PreProcessingImpl() {
  return true;
}

void VotincevDRadixMergeSortOMP::LocalRadixSort(uint32_t *begin, uint32_t *end) {
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

void VotincevDRadixMergeSortOMP::Merge(uint32_t *src, uint32_t *dst, int32_t left, int32_t mid, int32_t right) {
  int32_t i = left, j = mid, k = left;
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

bool VotincevDRadixMergeSortOMP::RunImpl() {
  const auto &input = GetInput();
  int32_t n = static_cast<int32_t>(input.size());
  if (n == 0) {
    return false;
  }

  //  uint32_t для работы с полным диапазоном int32_t без переполнений
  std::vector<uint32_t> working_array(n);
  int32_t min_val = input[0];

#pragma omp parallel for reduction(min : min_val)
  for (int32_t i = 0; i < n; ++i) {
    if (input[i] < min_val) {
      min_val = input[i];
    }
  }

#pragma omp parallel for
  for (int32_t i = 0; i < n; ++i) {
    working_array[i] = static_cast<uint32_t>(input[i]) - static_cast<uint32_t>(min_val);
  }

  std::vector<uint32_t> temp_buffer(n);

#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    int n_threads = omp_get_num_threads();

    // равномерное распределение нагрузки
    int32_t items = n / n_threads;
    int32_t rem = n % n_threads;
    int32_t l = tid * items + std::min(tid, rem);
    int32_t r = l + items + (tid < rem ? 1 : 0);

    if (l < r) {
      LocalRadixSort(working_array.data() + l, working_array.data() + r);
    }

    // слияние блоков
    for (int32_t step = 1; step < n_threads; step *= 2) {
#pragma omp barrier
      if (tid % (2 * step) == 0 && tid + step < n_threads) {
        int32_t m = (tid + step) * items + std::min(tid + step, rem);
        int32_t next_r = (tid + 2 * step) * items + std::min(tid + 2 * step, rem);
        if (next_r > n) {
          next_r = n;
        }

        Merge(working_array.data(), temp_buffer.data(), l, m, next_r);
        std::copy(temp_buffer.data() + l, temp_buffer.data() + next_r, working_array.data() + l);
      }
    }
  }

  std::vector<int32_t> result(n);
#pragma omp parallel for
  for (int32_t i = 0; i < n; ++i) {
    result[i] = static_cast<int32_t>(working_array[i] + static_cast<uint32_t>(min_val));
  }

  GetOutput() = std::move(result);
  return true;
}

bool VotincevDRadixMergeSortOMP::PostProcessingImpl() {
  return true;
}

}  // namespace votincev_d_radixmerge_sort
