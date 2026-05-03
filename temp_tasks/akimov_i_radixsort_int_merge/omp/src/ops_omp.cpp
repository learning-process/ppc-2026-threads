#include "akimov_i_radixsort_int_merge/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <vector>

#include "akimov_i_radixsort_int_merge/common/include/common.hpp"
#include "util/include/util.hpp"

namespace akimov_i_radixsort_int_merge {

namespace {
void CountingSortStep(std::vector<int>::iterator in_begin, std::vector<int>::iterator in_end,
                      std::vector<int>::iterator out_begin, size_t byte_index) {
  size_t shift = byte_index * 8;
  std::array<size_t, 256> count = {0};

  for (auto it = in_begin; it != in_end; ++it) {
    auto raw_val = static_cast<unsigned int>(*it);
    unsigned int byte_val = (raw_val >> shift) & 0xFF;

    if (byte_index == sizeof(int) - 1) {
      byte_val ^= 128;
    }
    count.at(byte_val)++;
  }

  std::array<size_t, 256> prefix{};
  prefix[0] = 0;
  for (int i = 1; i < 256; ++i) {
    prefix.at(i) = prefix.at(i - 1) + count.at(i - 1);
  }

  for (auto it = in_begin; it != in_end; ++it) {
    auto raw_val = static_cast<unsigned int>(*it);
    unsigned int byte_val = (raw_val >> shift) & 0xFF;

    if (byte_index == sizeof(int) - 1) {
      byte_val ^= 128;
    }

    *(out_begin + static_cast<int64_t>(prefix.at(byte_val))) = *it;
    prefix.at(byte_val)++;
  }
}

void RadixSortLocal(std::vector<int>::iterator begin, std::vector<int>::iterator end) {
  size_t n = std::distance(begin, end);
  if (n < 2) {
    return;
  }

  std::vector<int> temp(n);

  for (size_t i = 0; i < sizeof(int); ++i) {
    if (i % 2 == 0) {
      CountingSortStep(begin, end, temp.begin(), i);
    } else {
      CountingSortStep(temp.begin(), temp.end(), begin, i);
    }
  }
}
}  // namespace

AkimovIRadixSortIntMergeOMP::AkimovIRadixSortIntMergeOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType();
}

bool AkimovIRadixSortIntMergeOMP::ValidationImpl() {
  return !GetInput().empty();
}

bool AkimovIRadixSortIntMergeOMP::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool AkimovIRadixSortIntMergeOMP::RunImpl() {
  auto &arr = GetOutput();
  int n = static_cast<int>(arr.size());
  if (n == 0) {
    return true;
  }

  int num_threads = ppc::util::GetNumThreads();
  // Для маленьких массивов используем один поток
  if (n < num_threads * 100) {
    num_threads = 1;
  }

  if (num_threads == 1) {
    RadixSortLocal(arr.begin(), arr.end());
    return true;
  }

  // Разбиваем массив на блоки
  std::vector<int> offsets(num_threads + 1);
  int chunk = n / num_threads;
  int rem = n % num_threads;
  int curr = 0;
  for (int i = 0; i < num_threads; ++i) {
    offsets[i] = curr;
    curr += chunk + (i < rem ? 1 : 0);
  }
  offsets[num_threads] = n;

  // Параллельная сортировка блоков
#pragma omp parallel num_threads(num_threads) default(none) shared(offsets, arr)
  {
    int tid = omp_get_thread_num();
    auto begin = arr.begin() + offsets[tid];
    auto end = arr.begin() + offsets[tid + 1];
    RadixSortLocal(begin, end);
  }

  // Слияние блоков (параллельное на каждом шаге)
  for (int step = 1; step < num_threads; step *= 2) {
#pragma omp parallel for num_threads(num_threads) default(none) shared(step, num_threads, offsets, arr)
    for (int i = 0; i < num_threads; i += 2 * step) {
      if (i + step < num_threads) {
        auto begin = arr.begin() + offsets[i];
        auto middle = arr.begin() + offsets[i + step];
        int end_idx = std::min(i + (2 * step), num_threads);
        auto end = arr.begin() + offsets[end_idx];
        std::inplace_merge(begin, middle, end);
      }
    }
  }

  return true;
}

bool AkimovIRadixSortIntMergeOMP::PostProcessingImpl() {
  return true;
}

}  // namespace akimov_i_radixsort_int_merge
