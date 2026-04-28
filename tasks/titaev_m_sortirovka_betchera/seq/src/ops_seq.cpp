#include "titaev_m_sortirovka_betchera/seq/include/ops_seq.hpp"

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

TitaevSortirovkaBetcheraSEQ::TitaevSortirovkaBetcheraSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool TitaevSortirovkaBetcheraSEQ::ValidationImpl() {
  return !GetInput().empty();
}

bool TitaevSortirovkaBetcheraSEQ::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

void TitaevSortirovkaBetcheraSEQ::ConvertToKeys(const InType &input, std::vector<uint64_t> &keys) {
  const size_t n = input.size();
  for (size_t i = 0; i < n; i++) {
    keys[i] = DoubleToOrderedUint(input[i]);
  }
}

void TitaevSortirovkaBetcheraSEQ::RadixSort(std::vector<uint64_t> &keys) {
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

void TitaevSortirovkaBetcheraSEQ::ConvertFromKeys(const std::vector<uint64_t> &keys, OutType &output) {
  const size_t n = keys.size();
  output.resize(n);
  for (size_t i = 0; i < n; i++) {
    output[i] = OrderedUintToDouble(keys[i]);
  }
}

void TitaevSortirovkaBetcheraSEQ::BatcherStep(OutType &result, size_t n, size_t step, size_t stage) {
  for (size_t i = 0; i < n; i++) {
    size_t j = i ^ stage;
    if (j <= i || j >= n) {
      continue;
    }

    const bool ascending = (i & step) == 0;
    const bool need_swap = ascending ? result[i] > result[j] : result[i] < result[j];
    if (need_swap) {
      std::swap(result[i], result[j]);
    }
  }
}

void TitaevSortirovkaBetcheraSEQ::BatcherSort() {
  auto &result = GetOutput();
  const size_t n = result.size();

  for (size_t step = 1; step < n; step <<= 1) {
    for (size_t stage = step; stage > 0; stage >>= 1) {
      BatcherStep(result, n, step, stage);
    }
  }
}

bool TitaevSortirovkaBetcheraSEQ::RunImpl() {
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

bool TitaevSortirovkaBetcheraSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace titaev_m_sortirovka_betchera
