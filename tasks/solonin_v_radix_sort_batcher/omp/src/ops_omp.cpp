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
  const size_t kBase = 256;
  std::vector<int> freq(kBase, 0);
  std::vector<int> out(data.size());
  bool last = (pos == sizeof(int) - 1ULL);

  for (int v : data) {
    int bv = (v >> (pos * 8ULL)) & 0xFF;
    if (last) bv ^= 0x80;
    freq[bv]++;
  }
  for (size_t i = 1; i < kBase; ++i) {
    freq[i] += freq[i - 1];
  }
  for (int i = static_cast<int>(data.size()) - 1; i >= 0; --i) {
    int bv = (data[i] >> (pos * 8ULL)) & 0xFF;
    if (last) bv ^= 0x80;
    out[--freq[bv]] = data[i];
  }
  data = out;
}

void RadixSortBatcherOMP::SortChunk(std::vector<int> &data, int left, int right) {
  std::vector<int> chunk(data.begin() + left, data.begin() + right);
  for (size_t p = 0; p < sizeof(int); ++p) {
    SortByDigit(chunk, p);
  }
  std::copy(chunk.begin(), chunk.end(), data.begin() + left);
}

void RadixSortBatcherOMP::CompareSwapRange(std::vector<int> &data, int offset, int step, int half) {
  int lim = std::min(step, static_cast<int>(data.size()) - offset - step);
  for (int i = 0; i < lim; ++i) {
    int a = offset + i;
    int b = offset + i + step;
    if ((a & half) == (b & half) && data[a] > data[b]) {
      std::swap(data[a], data[b]);
    }
  }
}

void RadixSortBatcherOMP::BatcherNetwork(std::vector<int> &data, int p, int nthreads) {
  int n = static_cast<int>(data.size());
  for (int pv = p; pv < n; pv <<= 1) {
    int half = pv << 1;
    for (int step = pv; step > 0; step >>= 1) {
#pragma omp parallel for schedule(static) default(none) shared(n, pv, half, step, data) num_threads(nthreads)
      for (int j = step % pv; j < n - step; j += 2 * step) {
        CompareSwapRange(data, j, step, half);
      }
    }
  }
}

bool RadixSortBatcherOMP::ValidationImpl() {
  return !GetInput().empty();
}

bool RadixSortBatcherOMP::PreProcessingImpl() {
  return true;
}

bool RadixSortBatcherOMP::RunImpl() {
  GetOutput() = GetInput();
  int n = static_cast<int>(GetOutput().size());
  if (n <= 1) {
    return true;
  }

  int nthreads = ppc::util::GetNumThreads();
  int chunk = n / nthreads;

#pragma omp parallel default(none) shared(nthreads, chunk, n) num_threads(nthreads)
  {
    int tid = omp_get_thread_num();
    int lo = tid * chunk;
    int hi = (tid == nthreads - 1) ? n : lo + chunk;
    if (lo < hi) {
      SortChunk(GetOutput(), lo, hi);
    }
  }

  int start = 1;
  while (start < chunk) start <<= 1;
  BatcherNetwork(GetOutput(), start, nthreads);

  return true;
}

bool RadixSortBatcherOMP::PostProcessingImpl() {
  return true;
}

}  // namespace solonin_v_radix_sort_batcher
