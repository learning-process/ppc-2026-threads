#include "zenin_a_radix_sort_double_batcher_merge_omp/omp/include/ops_omp.hpp"

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <utility>
#include <vector>

#include "util/include/util.hpp"
#include "zenin_a_radix_sort_double_batcher_merge_omp/common/include/common.hpp"

namespace zenin_a_radix_sort_double_batcher_merge_omp {

ZeninARadixSortDoubleBatcherMergeOMP::ZeninARadixSortDoubleBatcherMergeOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool ZeninARadixSortDoubleBatcherMergeOMP::ValidationImpl() {
  return true;
}

bool ZeninARadixSortDoubleBatcherMergeOMP::PreProcessingImpl() {
  return true;
}

void ZeninARadixSortDoubleBatcherMergeOMP::BlocksComparing(std::vector<double> &arr, size_t i, size_t step) {
  for (size_t k = 0; k < step; ++k) {
    if (arr[i + k] > arr[i + k + step]) {
      std::swap(arr[i + k], arr[i + k + step]);
    }
  }
}

void ZeninARadixSortDoubleBatcherMergeOMP::BatcherOddEvenMerge(std::vector<double> &arr, size_t n) {
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

uint64_t ZeninARadixSortDoubleBatcherMergeOMP::PackDouble(double v) noexcept {
  uint64_t bits = 0ULL;
  std::memcpy(&bits, &v, sizeof(bits));
  if ((bits & (1ULL << 63)) != 0ULL) {
    bits = ~bits;
  } else {
    bits ^= (1ULL << 63);
  }
  return bits;
}

double ZeninARadixSortDoubleBatcherMergeOMP::UnpackDouble(uint64_t k) noexcept {
  if ((k & (1ULL << 63)) != 0ULL) {
    k ^= (1ULL << 63);
  } else {
    k = ~k;
  }
  double v = 0.0;
  std::memcpy(&v, &k, sizeof(v));
  return v;
}

void ZeninARadixSortDoubleBatcherMergeOMP::LSDRadixSort(std::vector<double> &array) {
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

bool ZeninARadixSortDoubleBatcherMergeOMP::RunImpl() {
  auto data = GetInput();
  if (data.empty()) {
    GetOutput() = data;
    return true;
  }
  size_t original_size = data.size();

  size_t pow2 = 1;
  while (pow2 < original_size) {
    pow2 *= 2;
  }
  data.resize(pow2, std::numeric_limits<double>::max());

  int num_threads = ppc::util::GetNumThreads();
  if (num_threads <= 0) {
    num_threads = 1;
  }

  size_t num_chunks = 1;
  while (num_chunks * 2 <= static_cast<size_t>(num_threads) && num_chunks * 2 <= pow2) {
    num_chunks *= 2;
  }

  size_t chunk_size = pow2 / num_chunks;
  int num_chunks_int = static_cast<int>(num_chunks);
  double *raw_data = data.data();

#pragma omp parallel for num_threads(num_threads) default(none) shared(chunk_size, raw_data, num_chunks_int)
  for (int i = 0; i < num_chunks_int; ++i) {
    size_t start = static_cast<size_t>(i) * chunk_size;
    std::vector<double> chunk(raw_data + start, raw_data + start + chunk_size);
    LSDRadixSort(chunk);
    for (size_t j = 0; j < chunk_size; ++j) {
      raw_data[start + j] = chunk[j];
    }
  }

  for (size_t size = chunk_size; size < pow2; size *= 2) {
    int merges_count = static_cast<int>(pow2 / (size * 2));
#pragma omp parallel for num_threads(num_threads) default(none) shared(size, raw_data, merges_count)
    for (int i = 0; i < merges_count; ++i) {
      size_t lo = static_cast<size_t>(i) * (2 * size);
      std::vector<double> block(raw_data + lo, raw_data + lo + (2 * size));
      BatcherOddEvenMerge(block, (2 * size));
      for (size_t j = 0; j < (2 * size); ++j) {
        raw_data[lo + j] = block[j];
      }
    }
  }

  data.resize(original_size);
  GetOutput() = data;
  return true;
}

bool ZeninARadixSortDoubleBatcherMergeOMP::PostProcessingImpl() {
  return true;
}

}  // namespace zenin_a_radix_sort_double_batcher_merge_omp
