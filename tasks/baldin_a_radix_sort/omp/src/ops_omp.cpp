#include "baldin_a_radix_sort/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <vector>

#include "baldin_a_radix_sort/common/include/common.hpp"
#include "util/include/util.hpp"

namespace baldin_a_radix_sort {

BaldinARadixSortOMP::BaldinARadixSortOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool BaldinARadixSortOMP::ValidationImpl() {
  return true;
}

bool BaldinARadixSortOMP::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

namespace {

void CountingSortStep(std::vector<int>::iterator in_begin, std::vector<int>::iterator in_end,
                      std::vector<int>::iterator out_begin, size_t byte_index) {
  size_t shift = byte_index * 8;
  std::array<size_t, 256> count = {0};

  for (auto it = in_begin; it != in_end; it++) {
    auto raw_val = static_cast<unsigned int>(*it);
    unsigned int byte_val = (raw_val >> shift) & 0xFF;

    if (byte_index == sizeof(int) - 1) {
      byte_val ^= 128;
    }
    count.at(byte_val)++;
  }

  std::array<size_t, 256> prefix{};
  prefix[0] = 0;
  for (int i = 1; i < 256; i++) {
    prefix.at(i) = prefix.at(i - 1) + count.at(i - 1);
  }

  for (auto it = in_begin; it != in_end; it++) {
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

  for (size_t i = 0; i < sizeof(int); i++) {
    size_t shift = i;

    if (i % 2 == 0) {
      CountingSortStep(begin, end, temp.begin(), shift);
    } else {
      CountingSortStep(temp.begin(), temp.end(), begin, shift);
    }
  }
}

}  // namespace

bool BaldinARadixSortOMP::RunImpl() {
  auto &out = GetOutput();
  int n = static_cast<int>(out.size());
  if (n == 0) {
    return true;
  }

  int num_threads = ppc::util::GetNumThreads();

  if (n < num_threads * 100) {
    num_threads = 1;
  }

  if (num_threads == 1) {
    RadixSortLocal(out.begin(), out.end());
    return true;
  }

  std::vector<int> offsets(num_threads + 1);
  int chunk = n / num_threads;
  int rem = n % num_threads;
  int curr = 0;
  for (int i = 0; i < num_threads; i++) {
    offsets[i] = curr;
    curr += chunk + (i < rem ? 1 : 0);
  }
  offsets[num_threads] = n;

#pragma omp parallel num_threads(num_threads) default(none) shared(num_threads, offsets, out)
  {
    int tid = omp_get_thread_num();
    auto begin = out.begin() + offsets[tid];
    auto end = out.begin() + offsets[tid + 1];

    RadixSortLocal(begin, end);
  }

  for (int step = 1; step < num_threads; step *= 2) {
#pragma omp parallel for num_threads(num_threads) default(none) shared(step, num_threads, offsets, out)
    for (int i = 0; i < num_threads; i += 2 * step) {
      if (i + step < num_threads) {
        auto begin = out.begin() + offsets[i];
        auto middle = out.begin() + offsets[i + step];

        int end_idx = std::min(i + (2 * step), num_threads);
        auto end = out.begin() + offsets[end_idx];

        std::inplace_merge(begin, middle, end);
      }
    }
  }

  return true;
}

bool BaldinARadixSortOMP::PostProcessingImpl() {
  return true;
}

}  // namespace baldin_a_radix_sort
