#include "zenin_a_radix_sort_double_batcher_merge/all/include/ops_all.hpp"

#include <mpi.h>
#include <tbb/parallel_invoke.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

#include "zenin_a_radix_sort_double_batcher_merge/common/include/common.hpp"

namespace zenin_a_radix_sort_double_batcher_merge {

ZeninARadixSortDoubleBatcherMergeALL::ZeninARadixSortDoubleBatcherMergeALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool ZeninARadixSortDoubleBatcherMergeALL::ValidationImpl() {
  return true;
}
bool ZeninARadixSortDoubleBatcherMergeALL::PreProcessingImpl() {
  local_data_ = GetInput();
  return true;
}

bool ZeninARadixSortDoubleBatcherMergeALL::PostProcessingImpl() {
  GetOutput() = local_data_;
  return true;
}
uint64_t ZeninARadixSortDoubleBatcherMergeALL::PackDouble(double v) noexcept {
  uint64_t bits = 0ULL;
  std::memcpy(&bits, &v, sizeof(bits));
  if ((bits & (1ULL << 63)) != 0ULL) {
    bits = ~bits;
  } else {
    bits ^= (1ULL << 63);
  }
  return bits;
}

double ZeninARadixSortDoubleBatcherMergeALL::UnpackDouble(uint64_t k) noexcept {
  if ((k & (1ULL << 63)) != 0ULL) {
    k ^= (1ULL << 63);
  } else {
    k = ~k;
  }
  double v = 0.0;
  std::memcpy(&v, &k, sizeof(v));
  return v;
}

void ZeninARadixSortDoubleBatcherMergeALL::LSDRadixSort(std::vector<double> &array) {
  const std::size_t n = array.size();
  if (n <= 1U) {
    return;
  }

  constexpr int kBits = 8;
  constexpr int kBuckets = 1 << kBits;
  constexpr int kPasses = static_cast<int>((sizeof(uint64_t) * 8) / kBits);

  std::vector<uint64_t> keys(n);
  for (std::size_t i = 0; i < n; ++i) {
    keys[i] = PackDouble(array[i]);
  }

  std::vector<uint64_t> tmp_keys(n);
  std::vector<double> tmp_vals(n);

  for (int pass = 0; pass < kPasses; ++pass) {
    int shift = pass * kBits;
    std::vector<std::size_t> cnt(kBuckets + 1, 0U);

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

void ZeninARadixSortDoubleBatcherMergeALL::BlocksComparing(std::vector<double> &arr, size_t i, size_t step) {
  for (size_t k = 0; k < step; ++k) {
    if (arr[i + k] > arr[i + k + step]) {
      std::swap(arr[i + k], arr[i + k + step]);
    }
  }
}

void ZeninARadixSortDoubleBatcherMergeALL::BatcherOddEvenMerge(std::vector<double> &arr, size_t n) {
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

void ZeninARadixSortDoubleBatcherMergeALL::SortChunk(std::vector<double> &chunk, size_t chunk_size) {
  size_t half = chunk_size / 2;
  if (half > 0) {
    std::vector<double> left(chunk.begin(), chunk.begin() + static_cast<std::ptrdiff_t>(half));
    std::vector<double> right(chunk.begin() + static_cast<std::ptrdiff_t>(half), chunk.end());
    tbb::parallel_invoke([&]() { LSDRadixSort(left); }, [&]() { LSDRadixSort(right); });
    std::ranges::copy(left, chunk.begin());
    std::ranges::copy(right, chunk.begin() + static_cast<std::ptrdiff_t>(half));
    BatcherOddEvenMerge(chunk, chunk_size);
  } else {
    LSDRadixSort(chunk);
  }
}

void ZeninARadixSortDoubleBatcherMergeALL::FinalMerge(size_t chunk_size, size_t pow2, size_t original_size) {
  for (size_t size = chunk_size; size < pow2; size *= 2) {
    size_t merges_count = pow2 / (size * 2);
    for (size_t i = 0; i < merges_count; ++i) {
      size_t lo = i * 2 * size;
      std::vector<double> block(local_data_.begin() + static_cast<std::ptrdiff_t>(lo),
                                local_data_.begin() + static_cast<std::ptrdiff_t>(lo + (2 * size)));
      BatcherOddEvenMerge(block, (2 * size));
      std::ranges::copy(block, local_data_.begin() + static_cast<std::ptrdiff_t>(lo));
    }
  }
  local_data_.resize(original_size);
}

bool ZeninARadixSortDoubleBatcherMergeALL::RunImpl() {
  int rank = 0;
  int num_procs = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  size_t original_size = 0;
  if (rank == 0) {
    original_size = local_data_.size();
  }
  MPI_Bcast(&original_size, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

  if (original_size == 0) {
    return true;
  }

  bool procs_is_pow2 = (num_procs > 1) && ((num_procs & (num_procs - 1)) == 0);

  if (!procs_is_pow2) {
    if (rank == 0) {
      LSDRadixSort(local_data_);
    }
    if (rank != 0) {
      local_data_.resize(original_size);
    }
    MPI_Bcast(local_data_.data(), static_cast<int>(original_size), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    return true;
  }

  size_t pow2 = 1;
  while (pow2 < original_size) {
    pow2 <<= 1;
  }
  while (pow2 % static_cast<size_t>(num_procs) != 0) {
    pow2 <<= 1;
  }
  MPI_Bcast(&pow2, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    local_data_.resize(pow2, std::numeric_limits<double>::max());
  }

  size_t chunk_size = pow2 / static_cast<size_t>(num_procs);
  std::vector<double> chunk(chunk_size);

  MPI_Scatter(local_data_.data(), static_cast<int>(chunk_size), MPI_DOUBLE, chunk.data(), static_cast<int>(chunk_size),
              MPI_DOUBLE, 0, MPI_COMM_WORLD);

  SortChunk(chunk, chunk_size);

  MPI_Gather(chunk.data(), static_cast<int>(chunk_size), MPI_DOUBLE, local_data_.data(), static_cast<int>(chunk_size),
             MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    FinalMerge(chunk_size, pow2, original_size);
  }

  if (rank != 0) {
    local_data_.resize(original_size);
  }
  MPI_Bcast(local_data_.data(), static_cast<int>(original_size), MPI_DOUBLE, 0, MPI_COMM_WORLD);

  return true;
}

}  // namespace zenin_a_radix_sort_double_batcher_merge
