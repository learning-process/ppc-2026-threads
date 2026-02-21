#include "zenin_a_radix_sort_double_batcher_merge_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

// #include "util/include/util.hpp"
#include "zenin_a_radix_sort_double_batcher_merge_seq/common/include/common.hpp"

namespace zenin_a_radix_sort_double_batcher_merge_seq {

ZeninARadixSortDoubleBatcherMergeSeqseq::ZeninARadixSortDoubleBatcherMergeSeqseq(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool ZeninARadixSortDoubleBatcherMergeSeqseq::ValidationImpl() {
  return true;
}

bool ZeninARadixSortDoubleBatcherMergeSeqseq::PreProcessingImpl() {
  return true;
}

void ZeninARadixSortDoubleBatcherMergeSeqseq::BlocksComparing(std::vector<double> &arr, size_t i, size_t step) {
  for (size_t k = 0; k < step; ++k) {
    if (arr[i + k] > arr[i + k + step]) {
      std::swap(arr[i + k], arr[i + k + step]);
    }
  }
}

void ZeninARadixSortDoubleBatcherMergeSeqseq::BatcherOddEvenMerge(std::vector<double> &arr, size_t n) {
  if (n <= 1) {
    return;
  }

  size_t step = n / 2;
  BlocksComparing(arr, 0, step);

  step /= 2;
  for (; step > 0; step /= 2) {
    for (size_t i = step; i < n - step; i += step * 2) {
      BlocksComparing(arr, i, step);
    }
  }
}

void ZeninARadixSortDoubleBatcherMergeSeqseq::BatcherMergeSort(std::vector<double> &array) {
  int n = static_cast<int>(array.size());
  if (n <= 1) {
    return;
  }

  int padded_n = 1;
  while (padded_n < n) {
    padded_n <<= 1;
  }
  array.resize(padded_n, std::numeric_limits<double>::max());
  int mid = padded_n / 2;

  std::vector<double> left(array.begin(), array.begin() + mid);
  std::vector<double> right(array.begin() + mid, array.end());

  LSDRadixSort(left);
  LSDRadixSort(right);

  std::ranges::copy(left.begin(), left.end(), array.begin());
  std::ranges::copy(right.begin(), right.end(), array.begin() + mid);

  BatcherOddEvenMerge(array, static_cast<size_t>(padded_n));
  array.resize(n);
}

uint64_t ZeninARadixSortDoubleBatcherMergeSeqseq::PackDouble(double v) noexcept {
  uint64_t bits = 0ULL;
  std::memcpy(&bits, &v, sizeof(bits));
  if ((bits & (1ULL << 63)) != 0ULL) {
    bits = ~bits;
  } else {
    bits ^= (1ULL << 63);
  }
  return bits;
}

double ZeninARadixSortDoubleBatcherMergeSeqseq::UnpackDouble(uint64_t k) noexcept {
  if ((k & (1ULL << 63)) != 0ULL) {
    k ^= (1ULL << 63);
  } else {
    k = ~k;
  }
  double v = 0.0;
  std::memcpy(&v, &k, sizeof(v));
  return v;
}

void ZeninARadixSortDoubleBatcherMergeSeqseq::LSDRadixSort(std::vector<double> &array) {
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

bool ZeninARadixSortDoubleBatcherMergeSeqseq::RunImpl() {
  std::vector<double> data = GetInput();
  BatcherMergeSort(data);
  GetOutput() = data;
  return true;
}

bool ZeninARadixSortDoubleBatcherMergeSeqseq::PostProcessingImpl() {
  return true;
}

}  // namespace zenin_a_radix_sort_double_batcher_merge_seq
