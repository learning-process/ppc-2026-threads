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
  uint64_t x = 0;
  std::memcpy(&x, &value, sizeof(double));
  constexpr uint64_t kSignMask = (1ULL << 63);
  if ((x & kSignMask) != 0ULL) {
    x = ~x;
  } else {
    x ^= kSignMask;
  }
  return x;
}

double OrderedUintToDouble(uint64_t x) {
  constexpr uint64_t kSignMask = (1ULL << 63);
  if ((x & kSignMask) != 0ULL) {
    x ^= kSignMask;
  } else {
    x = ~x;
  }
  double result = 0.0;
  std::memcpy(&result, &x, sizeof(double));
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
  const size_t n = input.size();
#pragma omp parallel for
  for (long long i = 0; i < static_cast<long long>(n); i++) {
    keys[i] = DoubleToOrderedUint(input[i]);
  }
}

void TitaevSortirovkaBetcheraALL::RadixSort(std::vector<uint64_t> &keys) {
  const size_t n = keys.size();
  if (n <= 1) {
    return;
  }

  constexpr int kBits = 8;
  constexpr int kBuckets = 1 << kBits;
  constexpr int kPasses = 64 / kBits;

  std::vector<uint64_t> tmp(n);

  for (int pass = 0; pass < kPasses; pass++) {
    std::vector<size_t> count(kBuckets, 0);
    const int num_threads = omp_get_max_threads();
    std::vector<std::vector<size_t>> local_count(num_threads, std::vector<size_t>(kBuckets, 0));

#pragma omp parallel
    {
      const int tid = omp_get_thread_num();
      auto &lc = local_count[tid];
#pragma omp for
      for (long long i = 0; i < static_cast<long long>(n); i++) {
        size_t bucket = (keys[i] >> (pass * kBits)) & (kBuckets - 1);
        lc[bucket]++;
      }
    }

    for (int b = 0; b < kBuckets; b++) {
      for (int t = 0; t < num_threads; t++) {
        count[b] += local_count[t][b];
      }
    }
    for (int i = 1; i < kBuckets; i++) {
      count[i] += count[i - 1];
    }
    for (size_t i = n; i-- > 0;) {
      size_t bucket = (keys[i] >> (pass * kBits)) & (kBuckets - 1);
      tmp[--count[bucket]] = keys[i];
    }
    keys.swap(tmp);
  }
}

void TitaevSortirovkaBetcheraALL::ConvertFromKeys(const std::vector<uint64_t> &keys, OutType &output) {
  const size_t n = keys.size();
  output.resize(n);
#pragma omp parallel for
  for (long long i = 0; i < static_cast<long long>(n); i++) {
    output[i] = OrderedUintToDouble(keys[i]);
  }
}

void TitaevSortirovkaBetcheraALL::BatcherSort() {
  auto &result = GetOutput();
  const size_t n = result.size();
  if (n < 2) {
    return;
  }

  for (size_t k = 2; k <= n; k <<= 1) {
    for (size_t j = k >> 1; j > 0; j >>= 1) {
#pragma omp parallel for
      for (long long ii = 0; ii < static_cast<long long>(n); ii++) {
        const size_t i = static_cast<size_t>(ii);
        const size_t l = i ^ j;
        if (l > i) {
          const bool ascending = ((i & k) == 0);
          const bool need_swap = ascending ? (result[i] > result[l]) : (result[i] < result[l]);
          if (need_swap) {
            std::swap(result[i], result[l]);
          }
        }
      }
    }
  }
}

bool TitaevSortirovkaBetcheraALL::RunImpl() {
  auto &input = GetInput();
  const size_t n = input.size();
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
