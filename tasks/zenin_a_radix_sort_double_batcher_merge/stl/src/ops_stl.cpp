#include "zenin_a_radix_sort_double_batcher_merge/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <thread>
#include <utility>
#include <vector>

#include "util/include/util.hpp"
#include "zenin_a_radix_sort_double_batcher_merge/common/include/common.hpp"

namespace zenin_a_radix_sort_double_batcher_merge {

ZeninARadixSortDoubleBatcherMergeSTL::ZeninARadixSortDoubleBatcherMergeSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool ZeninARadixSortDoubleBatcherMergeSTL::ValidationImpl() {
  return true;
}

bool ZeninARadixSortDoubleBatcherMergeSTL::PreProcessingImpl() {
  return true;
}

void ZeninARadixSortDoubleBatcherMergeSTL::BlocksComparing(std::vector<double> &arr, size_t i, size_t j) {
  if (arr[i] > arr[j]) {
    std::swap(arr[i], arr[j]);
  }
}

uint64_t ZeninARadixSortDoubleBatcherMergeSTL::PackDouble(double v) noexcept {
  uint64_t bits = 0ULL;
  std::memcpy(&bits, &v, sizeof(bits));
  if ((bits & (1ULL << 63)) != 0ULL) {
    bits = ~bits;
  } else {
    bits ^= (1ULL << 63);
  }
  return bits;
}

double ZeninARadixSortDoubleBatcherMergeSTL::UnpackDouble(uint64_t k) noexcept {
  if ((k & (1ULL << 63)) != 0ULL) {
    k ^= (1ULL << 63);
  } else {
    k = ~k;
  }
  double v = 0.0;
  std::memcpy(&v, &k, sizeof(v));
  return v;
}

void ZeninARadixSortDoubleBatcherMergeSTL::LSDRadixSort(std::vector<double> &array) {
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

void ZeninARadixSortDoubleBatcherMergeSTL::MergeFirstPass(std::vector<double> &arr, size_t po, size_t num_threads) {
  size_t chunk = (po + num_threads - 1) / num_threads;
  std::vector<std::thread> threads;
  for (size_t tid = 0; tid < num_threads; ++tid) {
    size_t lo = tid * chunk;
    size_t hi = std::min(lo + chunk, po);
    if (lo >= hi) {
      break;
    }
    threads.emplace_back([&arr, lo, hi, po]() {
      for (size_t i = lo; i < hi; ++i) {
        BlocksComparing(arr, i, i + po);
      }
    });
  }
  for (auto &th : threads) {
    th.join();
  }
}

void ZeninARadixSortDoubleBatcherMergeSTL::BatcherOddEvenMerge(std::vector<double> &arr, size_t n, size_t num_threads) {
  for (size_t po = n / 2; po > 0; po >>= 1) {
    if (po == n / 2) {
      MergeFirstPass(arr, po, num_threads);
    } else {
      for (size_t i = po; i < n - po; i += 2 * po) {
        for (size_t j = 0; j < po; ++j) {
          BlocksComparing(arr, i + j, i + j + po);
        }
      }
    }
  }
}

bool ZeninARadixSortDoubleBatcherMergeSTL::RunImpl() {
  auto data = GetInput();
  size_t original_size = data.size();

  if (original_size <= 1) {
    GetOutput() = data;
    return true;
  }

  int num_threads_int = ppc::util::GetNumThreads();
  if (num_threads_int <= 0) {
    num_threads_int = 1;
  }
  auto num_threads = static_cast<size_t>(num_threads_int);

  size_t pow2 = 1;
  while (pow2 < original_size) {
    pow2 <<= 1;
  }
  data.resize(pow2, std::numeric_limits<double>::max());

  size_t half = pow2 / 2;
  auto half_dist = static_cast<std::ptrdiff_t>(half);

  std::vector<double> left(data.begin(), data.begin() + half_dist);
  std::vector<double> right(data.begin() + half_dist, data.end());

  std::thread t([&]() { LSDRadixSort(left); });
  LSDRadixSort(right);
  t.join();

  std::ranges::copy(left, data.begin());
  std::ranges::copy(right, data.begin() + half_dist);

  BatcherOddEvenMerge(data, data.size(), num_threads);

  data.resize(original_size);
  GetOutput() = data;
  return true;
}

bool ZeninARadixSortDoubleBatcherMergeSTL::PostProcessingImpl() {
  return true;
}

}  // namespace zenin_a_radix_sort_double_batcher_merge
