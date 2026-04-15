#include "votincev_d_radixmerge_sort/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "votincev_d_radixmerge_sort/common/include/common.hpp"

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
      size_t digit = static_cast<size_t>((src[i] / exp) % 10);

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

void VotincevDRadixMergeSortOMP::Merge(const uint32_t *src, uint32_t *dst, int32_t left, int32_t mid, int32_t right) {
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

bool VotincevDRadixMergeSortOMP::RunImpl() {
  const auto &input = GetInput();
  auto n = static_cast<int32_t>(input.size());

  std::vector<uint32_t> working_array(static_cast<size_t>(n));
  int32_t min_val = input[0];

#pragma omp parallel for reduction(min : min_val) default(none) shared(n, input)
  for (int32_t i = 0; i < n; ++i) {
    min_val = std::min(input[i], min_val);
  }

#pragma omp parallel for default(none) shared(n, working_array, input, min_val)
  for (int32_t i = 0; i < n; ++i) {
    working_array[static_cast<size_t>(i)] = static_cast<uint32_t>(input[i]) - static_cast<uint32_t>(min_val);
  }

  std::vector<uint32_t> temp_buffer(static_cast<size_t>(n));

#pragma omp parallel default(none) shared(n, working_array, temp_buffer)
  {
    int tid = omp_get_thread_num();
    int n_threads = omp_get_num_threads();

    int32_t items = n / n_threads;
    int32_t rem = n % n_threads;
    int32_t l = (tid * items) + std::min(tid, rem);
    int32_t r = l + items + (tid < rem ? 1 : 0);

    if (l < r) {
      LocalRadixSort(working_array.data() + l, working_array.data() + r);
    }

    for (int32_t step = 1; step < n_threads; step *= 2) {
#pragma omp barrier
      if ((tid % (2 * step) == 0) && (tid + step < n_threads)) {
        int32_t m = ((tid + step) * items) + std::min(tid + step, rem);
        int32_t next_tid = std::min(tid + (2 * step), n_threads);
        int32_t next_r = (next_tid * items) + std::min(next_tid, rem);

        Merge(working_array.data(), temp_buffer.data(), l, m, next_r);
        std::copy(temp_buffer.data() + l, temp_buffer.data() + next_r, working_array.data() + l);
      }
    }
  }

  std::vector<int32_t> result(static_cast<size_t>(n));
#pragma omp parallel for default(none) shared(n, result, working_array, min_val)
  for (int32_t i = 0; i < n; ++i) {
    result[static_cast<size_t>(i)] =
        static_cast<int32_t>(working_array[static_cast<size_t>(i)] + static_cast<uint32_t>(min_val));
  }

  GetOutput() = std::move(result);
  return true;
}

bool VotincevDRadixMergeSortOMP::PostProcessingImpl() {
  return true;
}

}  // namespace votincev_d_radixmerge_sort
