#include "titaev_m_sortirovka_betchera/tbb/include/ops_tbb.hpp"

#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/parallel_for.h>

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

TitaevSortirovkaBetcheraTBB::TitaevSortirovkaBetcheraTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool TitaevSortirovkaBetcheraTBB::ValidationImpl() {
  return !GetInput().empty();
}

bool TitaevSortirovkaBetcheraTBB::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

void TitaevSortirovkaBetcheraTBB::ConvertToKeys(const InType &input, std::vector<uint64_t> &keys) {
  const size_t n = input.size();
  oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<size_t>(0, n), [&](const oneapi::tbb::blocked_range<size_t> &r) {
    for (size_t i = r.begin(); i < r.end(); i++) {
      keys[i] = DoubleToOrderedUint(input[i]);
    }
  });
}

void TitaevSortirovkaBetcheraTBB::RadixSort(std::vector<uint64_t> &keys) {
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
    for (size_t i = 0; i < n; i++) {
      size_t bucket = (keys[i] >> (pass * kBits)) & (kBuckets - 1);
      count[bucket]++;
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

void TitaevSortirovkaBetcheraTBB::ConvertFromKeys(const std::vector<uint64_t> &keys, OutType &output) {
  const size_t n = keys.size();
  output.resize(n);
  oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<size_t>(0, n), [&](const oneapi::tbb::blocked_range<size_t> &r) {
    for (size_t i = r.begin(); i < r.end(); i++) {
      output[i] = OrderedUintToDouble(keys[i]);
    }
  });
}

void TitaevSortirovkaBetcheraTBB::BatcherSort() {
  auto &result = GetOutput();
  const size_t n = result.size();
  if (n < 2) {
    return;
  }

  for (size_t k = 2; k <= n; k <<= 1) {
    for (size_t j = k >> 1; j > 0; j >>= 1) {
      oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<size_t>(0, n),
                                [&](const oneapi::tbb::blocked_range<size_t> &r) {
        for (size_t i = r.begin(); i < r.end(); i++) {
          const size_t l = i ^ j;
          if (l > i && l < n) {
            const bool ascending = ((i & k) == 0);
            const bool need_swap = ascending ? (result[i] > result[l]) : (result[i] < result[l]);
            if (need_swap) {
              std::swap(result[i], result[l]);
            }
          }
        }
      });
    }
  }
}

bool TitaevSortirovkaBetcheraTBB::RunImpl() {
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

bool TitaevSortirovkaBetcheraTBB::PostProcessingImpl() {
  return true;
}

}  // namespace titaev_m_sortirovka_betchera
