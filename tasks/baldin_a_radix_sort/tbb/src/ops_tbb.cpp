#include "baldin_a_radix_sort/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <vector>

#include "baldin_a_radix_sort/common/include/common.hpp"
#include "oneapi/tbb/parallel_for.h"
#include "util/include/util.hpp"

namespace baldin_a_radix_sort {

BaldinARadixSortTBB::BaldinARadixSortTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool BaldinARadixSortTBB::ValidationImpl() {
  return true;
}

bool BaldinARadixSortTBB::PreProcessingImpl() {
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

bool BaldinARadixSortTBB::RunImpl() {
  auto &out = GetOutput();
  int n = static_cast<int>(out.size());

  int num_chunks = ppc::util::GetNumThreads();

  std::vector<int> offsets(num_chunks + 1);
  int chunk_size = n / num_chunks;
  int rem = n % num_chunks;
  int curr = 0;
  for (int i = 0; i < num_chunks; i++) {
    offsets[i] = curr;
    curr += chunk_size + (i < rem ? 1 : 0);
  }
  offsets[num_chunks] = n;

  tbb::parallel_for(tbb::blocked_range<int>(0, num_chunks), [&](const tbb::blocked_range<int> &r) {
    for (int tid = r.begin(); tid != r.end(); tid++) {
      auto begin = out.begin() + offsets[tid];
      auto end = out.begin() + offsets[tid + 1];
      RadixSortLocal(begin, end);
    }
  });

  for (int step = 1; step < num_chunks; step *= 2) {
    int num_merges = (num_chunks + (2 * step) - 1) / (2 * step);
    tbb::parallel_for(tbb::blocked_range<int>(0, num_merges), [&](const tbb::blocked_range<int> &r) {
      for (int m_idx = r.begin(); m_idx != r.end(); m_idx++) {
        int i = m_idx * (2 * step);

        if (i + step < num_chunks) {
          auto begin = out.begin() + offsets[i];
          auto middle = out.begin() + offsets[i + step];
          int end_idx = std::min(i + (2 * step), num_chunks);
          auto end = out.begin() + offsets[end_idx];

          std::inplace_merge(begin, middle, end);
        }
      }
    });
  }

  return true;
}

bool BaldinARadixSortTBB::PostProcessingImpl() {
  return true;
}

}  // namespace baldin_a_radix_sort
