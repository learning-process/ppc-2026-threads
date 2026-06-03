#include "titaev_m_sortirovka_betchera/seq/include/ops_seq.hpp"

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
  const std::size_t n = input.size();
  for (std::size_t i = 0; i < n; i++) {
    keys[i] = DoubleToOrderedUint(input[i]);
  }
}

void TitaevSortirovkaBetcheraSEQ::RadixSort(std::vector<uint64_t> &keys) {
  const std::size_t n = keys.size();
  if (n <= 1) {
    return;
  }

  constexpr int kBits = 8;
  constexpr int kBuckets = 1 << kBits;
  constexpr int kPasses = 64 / kBits;

  std::vector<uint64_t> tmp(n);

  for (int pass = 0; pass < kPasses; pass++) {
    std::vector<std::size_t> count(kBuckets, 0);

    for (std::size_t i = 0; i < n; i++) {
      const std::size_t bucket = (keys[i] >> (pass * kBits)) & (kBuckets - 1);
      count[bucket]++;
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

void TitaevSortirovkaBetcheraSEQ::ConvertFromKeys(const std::vector<uint64_t> &keys, OutType &output) {
  const std::size_t n = keys.size();
  output.resize(n);
  for (std::size_t i = 0; i < n; i++) {
    output[i] = OrderedUintToDouble(keys[i]);
  }
}

void TitaevSortirovkaBetcheraSEQ::BatcherStage(OutType &result, std::size_t array_size, std::size_t block,
                                               std::size_t step) {
  for (std::size_t i = 0; i < array_size; i++) {
    const std::size_t partner = i ^ step;
    if (partner <= i) {
      continue;
    }
    const bool ascending = ((i & block) == 0);
    const bool need_swap = ascending ? (result[i] > result[partner]) : (result[i] < result[partner]);
    if (need_swap) {
      std::swap(result[i], result[partner]);
    }
  }
}

void TitaevSortirovkaBetcheraSEQ::BatcherSort() {
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

bool TitaevSortirovkaBetcheraSEQ::RunImpl() {
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

bool TitaevSortirovkaBetcheraSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace titaev_m_sortirovka_betchera
