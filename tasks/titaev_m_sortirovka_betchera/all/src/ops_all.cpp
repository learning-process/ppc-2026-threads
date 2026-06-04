#include "titaev_m_sortirovka_betchera/all/include/ops_all.hpp"

#include <omp.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <vector>

#include "titaev_m_sortirovka_betchera/common/include/common.hpp"

namespace titaev_m_sortirovka_betchera {

namespace {

uint64_t DoubleToOrderedUint(double value) {
  uint64_t bits = 0;
  std::memcpy(&bits, &value, sizeof(double));
  constexpr uint64_t kSignMask = (1ULL << 63);
  if ((bits & kSignMask) != 0ULL) {
    bits = ~bits;
  } else {
    bits ^= kSignMask;
  }
  return bits;
}

double OrderedUintToDouble(uint64_t bits) {
  constexpr uint64_t kSignMask = (1ULL << 63);
  if ((bits & kSignMask) != 0ULL) {
    bits ^= kSignMask;
  } else {
    bits = ~bits;
  }
  double result = 0.0;
  std::memcpy(&result, &bits, sizeof(double));
  return result;
}

}  // namespace

TitaevSortirovkaBetcheraALL::TitaevSortirovkaBetcheraALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool TitaevSortirovkaBetcheraALL::ValidationImpl() {
  return !GetInput().empty();
}

bool TitaevSortirovkaBetcheraALL::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

void TitaevSortirovkaBetcheraALL::ConvertToKeys(const InType &input, std::vector<uint64_t> &keys) {
  const auto n = static_cast<int64_t>(input.size());
#pragma omp parallel for default(none) shared(input, keys, n)
  for (int64_t i = 0; i < n; i++) {
    keys[i] = DoubleToOrderedUint(input[i]);
  }
}

void TitaevSortirovkaBetcheraALL::RadixSort(std::vector<uint64_t> &keys) {
  const std::size_t n = keys.size();
  if (n <= 1) {
    return;
  }

  constexpr int kBits = 8;
  constexpr int kBuckets = 1 << kBits;
  constexpr int kPasses = 64 / kBits;

  std::vector<uint64_t> tmp(n);
  const auto signed_n = static_cast<int64_t>(n);

  for (int pass = 0; pass < kPasses; pass++) {
    std::vector<std::size_t> count(kBuckets, 0);
    const int num_threads = omp_get_max_threads();
    std::vector<std::vector<std::size_t>> local_count(num_threads, std::vector<std::size_t>(kBuckets, 0));

#pragma omp parallel default(none) shared(keys, local_count, signed_n, pass)
    {
      const int tid = omp_get_thread_num();
      auto &lc = local_count[tid];
#pragma omp for
      for (int64_t i = 0; i < signed_n; i++) {
        const std::size_t bucket = (keys[i] >> (pass * kBits)) & (kBuckets - 1);
        lc[bucket]++;
      }
    }

    for (int bucket_idx = 0; bucket_idx < kBuckets; bucket_idx++) {
      for (int thr = 0; thr < num_threads; thr++) {
        count[bucket_idx] += local_count[thr][bucket_idx];
      }
    }

    for (int i = 1; i < kBuckets; i++) {
      count[i] += count[i - 1];
    }

    for (std::size_t i = n; i-- > 0;) {
      const std::size_t bucket = (keys[i] >> (pass * kBits)) & (kBuckets - 1);
      tmp[--count[bucket]] = keys[i];
    }

    keys.swap(tmp);
  }
}

void TitaevSortirovkaBetcheraALL::ConvertFromKeys(const std::vector<uint64_t> &keys, OutType &output) {
  const auto n = static_cast<int64_t>(keys.size());
  output.resize(keys.size());
#pragma omp parallel for default(none) shared(keys, output, n)
  for (int64_t i = 0; i < n; i++) {
    output[i] = OrderedUintToDouble(keys[i]);
  }
}

void TitaevSortirovkaBetcheraALL::BatcherStage(OutType &result, std::size_t array_size, std::size_t block,
                                               std::size_t step) {
  const auto signed_size = static_cast<int64_t>(array_size);
#pragma omp parallel for default(none) shared(result, signed_size, block, step)
  for (int64_t i = 0; i < signed_size; i++) {
    const auto idx = static_cast<std::size_t>(i);
    const std::size_t partner = idx ^ step;
    if (partner <= idx) {
      continue;
    }
    const bool ascending = ((idx & block) == 0);
    const bool need_swap = ascending ? (result[idx] > result[partner]) : (result[idx] < result[partner]);
    if (need_swap) {
      std::swap(result[idx], result[partner]);
    }
  }
}

void TitaevSortirovkaBetcheraALL::BatcherSort() {
  auto &result = GetOutput();
  const std::size_t n = result.size();
  if (n < 2) {
    return;
  }
  for (std::size_t block = 2; block <= n; block <<= 1) {
    for (std::size_t step = block >> 1; step > 0; step >>= 1) {
      BatcherStage(result, n, block, step);
    }
  }
}

bool TitaevSortirovkaBetcheraALL::RunImpl() {
  auto &input = GetInput();
  const std::size_t n = input.size();
  if (n <= 1) {
    return true;
  }
  std::vector<uint64_t> keys(n);
  ConvertToKeys(input, keys);
  RadixSort(keys);
  ConvertFromKeys(keys, GetOutput());
  if ((n & (n - 1)) == 0) {
    BatcherSort();
  }
  return true;
}

bool TitaevSortirovkaBetcheraALL::PostProcessingImpl() {
  return true;
}

}  // namespace titaev_m_sortirovka_betchera
