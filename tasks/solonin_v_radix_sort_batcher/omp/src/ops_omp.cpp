#include "solonin_v_radix_sort_batcher/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <vector>

#include "solonin_v_radix_sort_batcher/common/include/common.hpp"
#include "util/include/util.hpp"

namespace solonin_v_radix_sort_batcher {

RadixSortBatcherOMP::RadixSortBatcherOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

void RadixSortBatcherOMP::SortByDigit(std::vector<int> &data, size_t pos) {
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

void RadixSortBatcherOMP::SortChunk(std::vector<int> &data, int left, int right) {
  std::vector<int> chunk(data.begin() + left, data.begin() + right);
  for (size_t pos_idx = 0; pos_idx < sizeof(int); ++pos_idx) {
    SortByDigit(chunk, pos_idx);
  }
  std::ranges::copy(chunk, data.begin() + left);
}

bool RadixSortBatcherOMP::ValidationImpl() {
  return !GetInput().empty();
}

bool RadixSortBatcherOMP::PreProcessingImpl() {
  return true;
}

bool RadixSortBatcherOMP::RunImpl() {
  std::vector<int> tmp(GetInput().begin(), GetInput().end());
  int n = static_cast<int>(tmp.size());
  if (n <= 1) {
    GetOutput() = std::move(tmp);
    return true;
  }
  int nthreads = std::min(ppc::util::GetNumThreads(), n);
  int chunk = (n + nthreads - 1) / nthreads;

  std::vector<std::vector<int>> blocks;
  for (int i = 0; i < nthreads; ++i) {
    int lo = i * chunk;
    int hi = std::min(lo + chunk, n);
    if (lo >= n) {
      break;
    }
    blocks.emplace_back(tmp.begin() + lo, tmp.begin() + hi);
  }

#pragma omp parallel for schedule(static) default(none) shared(blocks) num_threads(nthreads)
  for (int i = 0; i < static_cast<int>(blocks.size()); ++i) {
    for (size_t pos_idx = 0; pos_idx < sizeof(int); ++pos_idx) {
      SortByDigit(blocks[i], pos_idx);
    }
  }

  while (blocks.size() > 1) {
    std::vector<std::vector<int>> merged;
    for (size_t i = 0; i + 1 < blocks.size(); i += 2) {
      std::vector<int> result;
      result.reserve(blocks[i].size() + blocks[i + 1].size());
      std::merge(blocks[i].begin(), blocks[i].end(), blocks[i + 1].begin(), blocks[i + 1].end(),
                 std::back_inserter(result));
      merged.push_back(std::move(result));
    }
    if (blocks.size() % 2 == 1) {
      merged.push_back(std::move(blocks.back()));
    }
    blocks = std::move(merged);
  }

  GetOutput() = std::move(blocks.front());
  return true;
}

bool RadixSortBatcherOMP::PostProcessingImpl() {
  return true;
}

}  // namespace solonin_v_radix_sort_batcher
