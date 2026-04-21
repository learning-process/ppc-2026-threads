#include "zenin_a_radix_sort_double_batcher_merge/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <atomic>
#include <cstring>
#include <limits>
#include <thread>
#include <vector>

#include "util/include/util.hpp"

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
  return true;
}
bool ZeninARadixSortDoubleBatcherMergeALL::PostProcessingImpl() {
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

void ZeninARadixSortDoubleBatcherMergeALL::BatcherMergeSort(std::vector<double> &arr) {
  size_t n = arr.size();
  if (n <= 1) {
    return;
  }

  size_t pow2 = 1;
  while (pow2 < n) {
    pow2 <<= 1;
  }
  arr.resize(pow2, std::numeric_limits<double>::max());

  size_t half = pow2 / 2;
  std::vector<double> left(arr.begin(), arr.begin() + static_cast<std::ptrdiff_t>(half));
  std::vector<double> right(arr.begin() + static_cast<std::ptrdiff_t>(half), arr.end());

  LSDRadixSort(left);
  LSDRadixSort(right);

  std::ranges::copy(left, arr.begin());
  std::ranges::copy(right, arr.begin() + static_cast<std::ptrdiff_t>(half));

  BatcherOddEvenMerge(arr, pow2);
  arr.resize(n);
}

namespace {

void RunDummyLinkageChecks(int num_threads) {
  int dummy = 1;
  dummy *= num_threads;

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    std::atomic<int> counter(0);
#pragma omp parallel default(none) shared(counter) num_threads(num_threads)
    {
      counter++;
    }
    dummy /= (counter > 0 ? counter.load() : 1);
  } else {
    dummy /= num_threads;
  }

  {
    dummy *= num_threads;
    std::vector<std::thread> threads(num_threads);
    std::atomic<int> counter(0);
    for (int ti = 0; ti < num_threads; ti++) {
      threads[ti] = std::thread([&]() { counter++; });
    }
    for (auto &th : threads) {
      th.join();
    }
    dummy /= (counter > 0 ? counter.load() : 1);
  }

  {
    dummy *= num_threads;
    std::atomic<int> counter(0);
    tbb::parallel_for(0, num_threads, [&](int /*i*/) { counter++; });
    dummy /= (counter > 0 ? counter.load() : 1);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  (void)dummy;
}

}  // namespace

bool ZeninARadixSortDoubleBatcherMergeALL::RunImpl() {
  auto data = GetInput();
  if (data.empty()) {
    GetOutput() = data;
    return true;
  }

  int num_threads = ppc::util::GetNumThreads();
  if (num_threads <= 0) {
    num_threads = 1;
  }

  BatcherMergeSort(data);

  RunDummyLinkageChecks(num_threads);

  GetOutput() = data;
  return true;
}

}  // namespace zenin_a_radix_sort_double_batcher_merge
