#include "solonin_v_radix_sort_batcher/tbb/include/ops_tbb.hpp"
#include "solonin_v_radix_sort_batcher/common/include/common.hpp"
#include "util/include/util.hpp"
#include <algorithm>
#include <cstddef>
#include <iterator>
#include <tbb/blocked_range.h>
#include <tbb/global_control.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_invoke.h>
#include <utility>
#include <vector>

namespace solonin_v_radix_sort_batcher {

RadixSortBatcherTBB::RadixSortBatcherTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

void RadixSortBatcherTBB::SortByDigit(std::vector<int> &data, size_t pos) {
  const size_t k_base = 256;
  std::vector<int> freq(k_base, 0);
  std::vector<int> out(data.size());
  bool last = (pos == sizeof(int) - 1ULL);

  for (int v : data) {
    int bv = (v >> (pos * 8ULL)) & 0xFF;
    if (last) {
      bv ^= 0x80;
    }
    freq[bv]++;
  }
  for (size_t i = 1; i < k_base; ++i) {
    freq[i] += freq[i - 1];
  }
  for (int i = static_cast<int>(data.size()) - 1; i >= 0; --i) {
    int bv = (data[i] >> (pos * 8ULL)) & 0xFF;
    if (last) {
      bv ^= 0x80;
    }
    out[--freq[bv]] = data[i];
  }
  data = out;
}

void RadixSortBatcherTBB::SortSegment(std::vector<int> &data, int lo, int hi) {
  std::vector<int> seg(data.begin() + lo, data.begin() + hi);
  for (size_t pos_idx = 0; pos_idx < sizeof(int); ++pos_idx) {
    SortByDigit(seg, pos_idx);
  }
  std::ranges::copy(seg, data.begin() + lo);
}

std::vector<int> RadixSortBatcherTBB::MergeBatcher(const std::vector<int> &left, const std::vector<int> &right) {
  if (left.empty()) {
    return right;
  }
  if (right.empty()) {
    return left;
  }

  std::vector<int> even_l;
  std::vector<int> odd_l;
  std::vector<int> even_r;
  std::vector<int> odd_r;

  for (size_t i = 0; i < left.size(); ++i) {
    if (i % 2 == 0) {
      even_l.push_back(left[i]);
    } else {
      odd_l.push_back(left[i]);
    }
  }
  for (size_t i = 0; i < right.size(); ++i) {
    if ((left.size() + i) % 2 == 0) {
      even_r.push_back(right[i]);
    } else {
      odd_r.push_back(right[i]);
    }
  }

  std::vector<int> merged_even;
  std::vector<int> merged_odd;

  tbb::parallel_invoke(
      [&]() { std::ranges::merge(even_l, even_r, std::back_inserter(merged_even)); },
      [&]() { std::ranges::merge(odd_l, odd_r, std::back_inserter(merged_odd)); });

  std::vector<int> result(left.size() + right.size());
  for (size_t i = 0; i < merged_even.size(); ++i) {
    result[i * 2] = merged_even[i];
  }
  for (size_t i = 0; i < merged_odd.size(); ++i) {
    result[(i * 2) + 1] = merged_odd[i];
  }

  for (size_t i = 1; i + 1 < result.size(); i += 2) {
    if (result[i] > result[i + 1]) {
      std::swap(result[i], result[i + 1]);
    }
  }

  return result;
}

bool RadixSortBatcherTBB::ValidationImpl() {
  return !GetInput().empty();
}

bool RadixSortBatcherTBB::PreProcessingImpl() {
  return true;
}

bool RadixSortBatcherTBB::RunImpl() {
  if (GetInput().size() <= 1) {
    GetOutput() = GetInput();
    return true;
  }

  int nthreads = std::max(1, ppc::util::GetNumThreads());
  const tbb::global_control ctrl(tbb::global_control::max_allowed_parallelism, nthreads);

  int n = static_cast<int>(GetInput().size());
  int chunk = (n + nthreads - 1) / nthreads;

  int actual_blocks = std::min(nthreads, (n + chunk - 1) / chunk);
  std::vector<std::vector<int>> blocks(actual_blocks);

  for (int thread_idx = 0; thread_idx < actual_blocks; ++thread_idx) {
    int lo = thread_idx * chunk;
    int hi = std::min(lo + chunk, n);
    blocks[thread_idx].assign(GetInput().begin() + lo, GetInput().begin() + hi);
  }

  tbb::parallel_for(tbb::blocked_range<int>(0, actual_blocks), [&](const tbb::blocked_range<int> &rng) {
    for (int thread_idx = rng.begin(); thread_idx < rng.end(); ++thread_idx) {
      for (size_t pos_idx = 0; pos_idx < sizeof(int); ++pos_idx) {
        SortByDigit(blocks[thread_idx], pos_idx);
      }
    }
  });

  while (blocks.size() > 1) {
    std::vector<std::vector<int>> next;
    for (size_t i = 0; i + 1 < blocks.size(); i += 2) {
      next.push_back(MergeBatcher(blocks[i], blocks[i + 1]));
    }
    if (blocks.size() % 2 == 1) {
      next.push_back(blocks.back());
    }
    blocks = std::move(next);
  }

  GetOutput() = blocks.front();
  return true;
}

bool RadixSortBatcherTBB::PostProcessingImpl() {
  return true;
}

}  // namespace solonin_v_radix_sort_batcher
