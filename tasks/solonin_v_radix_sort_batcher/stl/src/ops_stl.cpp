#include "solonin_v_radix_sort_batcher/stl/include/ops_stl.hpp"
#include <algorithm>
#include <cstddef>
#include <thread>
#include <vector>
#include "solonin_v_radix_sort_batcher/common/include/common.hpp"
#include "util/include/util.hpp"

namespace solonin_v_radix_sort_batcher {

RadixSortBatcherSTL::RadixSortBatcherSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

void RadixSortBatcherSTL::SortByDigit(std::vector<int> &data, size_t pos) {
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

void RadixSortBatcherSTL::SortBlock(std::vector<int> &data, int lo, int hi) {
  std::vector<int> seg(data.begin() + lo, data.begin() + hi);
  for (size_t p = 0; p < sizeof(int); ++p) {
    SortByDigit(seg, p);
  }
  std::copy(seg.begin(), seg.end(), data.begin() + lo);
}

std::vector<int> RadixSortBatcherSTL::OddEvenMerge(const std::vector<int> &a, const std::vector<int> &b) {
  if (a.empty()) return b;
  if (b.empty()) return a;

  std::vector<int> evens;
  std::vector<int> odds;

  for (size_t i = 0; i < a.size(); ++i) {
    if (i % 2 == 0) evens.push_back(a[i]);
    else odds.push_back(a[i]);
  }
  for (size_t i = 0; i < b.size(); ++i) {
    if ((a.size() + i) % 2 == 0) evens.push_back(b[i]);
    else odds.push_back(b[i]);
  }

  std::vector<int> se;
  std::vector<int> so;

  std::thread t_even([&]() { std::ranges::sort(evens); se = evens; });
  std::thread t_odd([&]() { std::ranges::sort(odds); so = odds; });
  t_even.join();
  t_odd.join();

  std::vector<int> result(a.size() + b.size());
  for (size_t i = 0; i < se.size(); ++i) result[i * 2] = se[i];
  for (size_t i = 0; i < so.size(); ++i) result[i * 2 + 1] = so[i];

  for (size_t i = 1; i + 1 < result.size(); i += 2) {
    if (result[i] > result[i + 1]) std::swap(result[i], result[i + 1]);
  }
  return result;
}

bool RadixSortBatcherSTL::ValidationImpl() {
  return !GetInput().empty();
}

bool RadixSortBatcherSTL::PreProcessingImpl() {
  return true;
}

bool RadixSortBatcherSTL::RunImpl() {
  if (GetInput().size() <= 1) {
    GetOutput() = GetInput();
    return true;
  }

  int nthreads = ppc::util::GetNumThreads();
  int n = static_cast<int>(GetInput().size());
  int chunk = (n + nthreads - 1) / nthreads;

  std::vector<std::vector<int>> blocks(nthreads);
  std::vector<std::thread> workers;

  for (int t = 0; t < nthreads; ++t) {
    int lo = t * chunk;
    int hi = std::min(lo + chunk, n);
    if (lo >= n) {
      blocks[t] = {};
      continue;
    }
    blocks[t].assign(GetInput().begin() + lo, GetInput().begin() + hi);
    workers.emplace_back([&blocks, t]() {
      for (size_t p = 0; p < sizeof(int); ++p) {
        SortByDigit(blocks[t], p);
      }
    });
  }
  for (auto &w : workers) w.join();

  while (blocks.size() > 1) {
    std::vector<std::vector<int>> merged;
    std::vector<std::thread> mt;
    merged.resize((blocks.size() + 1) / 2);

    for (size_t i = 0; i + 1 < blocks.size(); i += 2) {
      mt.emplace_back([&blocks, &merged, i]() {
        merged[i / 2] = OddEvenMerge(blocks[i], blocks[i + 1]);
      });
    }
    if (blocks.size() % 2 == 1) merged.back() = blocks.back();

    for (auto &m : mt) m.join();
    blocks = std::move(merged);
  }

  GetOutput() = blocks.front();
  return true;
}

bool RadixSortBatcherSTL::PostProcessingImpl() {
  return true;
}

}  // namespace solonin_v_radix_sort_batcher
