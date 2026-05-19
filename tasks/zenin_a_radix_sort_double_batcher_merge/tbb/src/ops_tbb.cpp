#include "zenin_a_radix_sort_double_batcher_merge/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <utility>
#include <vector>

#include "oneapi/tbb/parallel_for.h"
#include "zenin_a_radix_sort_double_batcher_merge/common/include/common.hpp"

namespace zenin_a_radix_sort_double_batcher_merge {

ZeninARadixSortDoubleBatcherMergeTBB::ZeninARadixSortDoubleBatcherMergeTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool ZeninARadixSortDoubleBatcherMergeTBB::ValidationImpl() {
  return true;
}

bool ZeninARadixSortDoubleBatcherMergeTBB::PreProcessingImpl() {
  return true;
}

void ZeninARadixSortDoubleBatcherMergeTBB::BlocksComparing(std::vector<double> &arr, size_t i, size_t j) {
  if (arr[i] > arr[j]) {
    std::swap(arr[i], arr[j]);
  }
}

uint64_t ZeninARadixSortDoubleBatcherMergeTBB::PackDouble(double v) noexcept {
  uint64_t bits = 0ULL;
  std::memcpy(&bits, &v, sizeof(bits));
  if ((bits & (1ULL << 63)) != 0ULL) {
    bits = ~bits;
  } else {
    bits ^= (1ULL << 63);
  }
  return bits;
}

double ZeninARadixSortDoubleBatcherMergeTBB::UnpackDouble(uint64_t k) noexcept {
  if ((k & (1ULL << 63)) != 0ULL) {
    k ^= (1ULL << 63);
  } else {
    k = ~k;
  }
  double v = 0.0;
  std::memcpy(&v, &k, sizeof(v));
  return v;
}

void ZeninARadixSortDoubleBatcherMergeTBB::LSDRadixSort(std::vector<double> &array) {
  const std::size_t n = array.size();
  if (n <= 1U) {
    return;
  }

  constexpr int kBits = 8;
  constexpr int kBuckets = 1 << kBits;
  constexpr int kPasses = static_cast<int>((sizeof(uint64_t) * 8) / kBits);

  std::vector<uint64_t> keys;
  keys.resize(n);
  for (std::size_t i = 0; i < n; ++i) {
    keys[i] = PackDouble(array[i]);
  }

  std::vector<uint64_t> tmp_keys;
  tmp_keys.resize(n);
  std::vector<double> tmp_vals;
  tmp_vals.resize(n);

  for (int pass = 0; pass < kPasses; ++pass) {
    int shift = pass * kBits;
    std::vector<std::size_t> cnt;
    cnt.assign(kBuckets + 1, 0U);

    for (std::size_t i = 0; i < n; ++i) {
      auto d = static_cast<std::size_t>((keys[i] >> shift) & (kBuckets - 1));
      ++cnt[d + 1];
    }
    for (int i = 0; i < kBuckets; ++i) {
      cnt[i + 1] += cnt[i];
    }

    for (std::size_t i = 0; i < n; ++i) {
      auto d = static_cast<std::size_t>((keys[i] >> shift) & (kBuckets - 1));
      std::size_t pos = cnt[d]++;
      tmp_keys[pos] = keys[i];
      tmp_vals[pos] = array[i];
    }

    keys.swap(tmp_keys);
    array.swap(tmp_vals);
  }

  for (std::size_t i = 0; i < n; ++i) {
    array[i] = UnpackDouble(keys[i]);
  }
}

void ZeninARadixSortDoubleBatcherMergeTBB::BatcherOddEvenMerge(std::vector<double> &arr, size_t n) {
  for (size_t po = n / 2; po > 0; po >>= 1) {
    if (po == n / 2) {
      tbb::parallel_for(tbb::blocked_range<size_t>(0, po), [&](const tbb::blocked_range<size_t> &r) {
        for (size_t i = r.begin(); i < r.end(); ++i) {
          BlocksComparing(arr, i, i + po);
        }
      });
    } else {
      for (size_t i = po; i < n - po; i += 2 * po) {
        for (size_t j = 0; j < po; ++j) {
          BlocksComparing(arr, i + j, i + j + po);
        }
      }
    }
  }
}

bool ZeninARadixSortDoubleBatcherMergeTBB::RunImpl() {
  auto data = GetInput();
  size_t original_size = data.size();

  if (original_size <= 1) {
    GetOutput() = data;
    return true;
  }

  size_t pow2 = 1;
  while (pow2 < original_size) {
    pow2 <<= 1;
  }
  data.resize(pow2, std::numeric_limits<double>::max());

  size_t half = pow2 / 2;
  auto half_dist = static_cast<std::ptrdiff_t>(half);

  std::vector<double> left(data.begin(), data.begin() + half_dist);
  std::vector<double> right(data.begin() + half_dist, data.end());

  tbb::parallel_invoke([&]() { LSDRadixSort(left); }, [&]() { LSDRadixSort(right); });

  std::ranges::copy(left, data.begin());
  std::ranges::copy(right, data.begin() + half_dist);

  BatcherOddEvenMerge(data, data.size());

  data.resize(original_size);
  GetOutput() = data;
  return true;
}

bool ZeninARadixSortDoubleBatcherMergeTBB::PostProcessingImpl() {
  return true;
}

}  // namespace zenin_a_radix_sort_double_batcher_merge
