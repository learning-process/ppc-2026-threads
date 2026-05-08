#include "gonozov_l_bitwise_sorting_double_Batcher_merge/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <execution>
#include <limits>
#include <numeric>
#include <thread>
#include <vector>

#include "gonozov_l_bitwise_sorting_double_Batcher_merge/common/include/common.hpp"

namespace gonozov_l_bitwise_sorting_double_batcher_merge {

GonozovLBitSortBatcherMergeSTL::GonozovLBitSortBatcherMergeSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool GonozovLBitSortBatcherMergeSTL::ValidationImpl() {
  return !GetInput().empty();  // проверка на то, что исходный массив непустой
}

bool GonozovLBitSortBatcherMergeSTL::PreProcessingImpl() {
  return true;
}

namespace {

constexpr size_t kRadix = 256;
constexpr size_t kCutoff = 1 << 18;

inline uint64_t ToKey(double d) {
  uint64_t x;
  std::memcpy(&x, &d, sizeof(double));
  return (x & (1ULL << 63)) ? ~x : (x | 0x8000000000000000ULL);
}

inline double ToDouble(uint64_t x) {
  x = (x & (1ULL << 63)) ? (x & ~0x8000000000000000ULL) : ~x;
  double d;
  std::memcpy(&d, &x, sizeof(double));
  return d;
}

size_t NextPow2(size_t n) {
  size_t p = 1;
  while (p < n) {
    p <<= 1;
  }
  return p;
}

thread_local std::vector<uint64_t> tl_a;
thread_local std::vector<uint64_t> tl_b;

void Radix(std::vector<double> &a, size_t l, size_t r) {
  size_t n = r - l;
  if (n <= 1) {
    return;
  }

  tl_a.resize(n);
  tl_b.resize(n);

  for (size_t i = 0; i < n; ++i) {
    tl_a[i] = ToKey(a[l + i]);
  }

  for (int p = 0; p < 8; ++p) {
    size_t cnt[kRadix] = {};
    int shift = p * 8;

    for (size_t i = 0; i < n; ++i) {
      cnt[(tl_a[i] >> shift) & 0xFF]++;
    }

    for (size_t i = 1; i < kRadix; ++i) {
      cnt[i] += cnt[i - 1];
    }

    for (int i = (int)n - 1; i >= 0; --i) {
      uint8_t b = (tl_a[i] >> shift) & 0xFF;
      tl_b[--cnt[b]] = tl_a[i];
    }

    tl_a.swap(tl_b);
  }

  for (size_t i = 0; i < n; ++i) {
    a[l + i] = ToDouble(tl_a[i]);
  }
}

inline void Merge(std::vector<double> &a, size_t i, size_t len) {
  size_t half = len / 2;
  size_t end = std::min(i + len, a.size());

  for (size_t s = half; s > 0; s >>= 1) {
    for (size_t j = i; j + s < end; ++j) {
      if (a[j] > a[j + s]) {
        std::swap(a[j], a[j + s]);
      }
    }
  }
}

void Batcher(std::vector<double> &a, size_t n) {
  for (size_t len = 2; len <= n; len <<= 1) {
    for (size_t i = 0; i < n; i += len) {
      Merge(a, i, len);
    }
  }
}

void HybridSortDouble(std::vector<double> &a) {
  if (a.size() <= 1) {
    return;
  }

  size_t n0 = a.size();
  size_t n = NextPow2(n0);

  a.resize(n, std::numeric_limits<double>::infinity());

  unsigned hw = std::thread::hardware_concurrency();
  if (hw == 0) {
    hw = 4;
  }

  size_t chunk = std::max(n / hw, kCutoff);

  std::vector<std::thread> threads;

  for (size_t t = 0; t < hw; ++t) {
    size_t l = t * chunk;
    if (l >= n) {
      break;
    }

    size_t r = std::min(l + chunk, n);

    threads.emplace_back([&, l, r] { Radix(a, l, r); });
  }

  for (auto &th : threads) {
    th.join();
  }

  Batcher(a, n);
  a.resize(n0);
}

}  // namespace

bool GonozovLBitSortBatcherMergeSTL::RunImpl() {
  std::vector<double> array = GetInput();
  HybridSortDouble(array);
  GetOutput() = array;
  return true;
}

bool GonozovLBitSortBatcherMergeSTL::PostProcessingImpl() {
  return true;
}

}  // namespace gonozov_l_bitwise_sorting_double_batcher_merge
