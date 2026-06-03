#include "titaev_m_sortirovka_betchera/omp/include/ops_omp.hpp"

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

void RadixPass(int pass, size_t n, const std::vector<uint64_t> &src, std::vector<uint64_t> &dest) {
  constexpr int kBits = 8;
  constexpr int kBuckets = 256;
  int num_threads = omp_get_max_threads();
  std::vector<std::vector<size_t>> local_counts(num_threads, std::vector<size_t>(kBuckets, 0));

#pragma omp parallel default(none) shared(src, local_counts, n, pass, kBits)
  {
    int tid = omp_get_thread_num();
#pragma omp for
    for (size_t i = 0; i < n; ++i) {
      size_t bucket = (src[i] >> (static_cast<size_t>(pass) * kBits)) & 255;
      local_counts[static_cast<size_t>(tid)][bucket]++;
    }
  }

  std::vector<size_t> common_counts(kBuckets, 0);
  for (const auto &l_count : local_counts) {
    for (size_t b_idx = 0; b_idx < kBuckets; ++b_idx) {
      common_counts[b_idx] += l_count[b_idx];
    }
  }

  std::vector<size_t> prefixes(kBuckets, 0);
  for (size_t b_idx = 1; b_idx < kBuckets; ++b_idx) {
    prefixes[b_idx] = prefixes[b_idx - 1] + common_counts[b_idx - 1];
  }

  std::vector<std::vector<size_t>> offsets(num_threads, std::vector<size_t>(kBuckets));
  for (size_t b_idx = 0; b_idx < kBuckets; ++b_idx) {
    size_t curr = prefixes[b_idx];
    for (int t_idx = 0; t_idx < num_threads; ++t_idx) {
      offsets[static_cast<size_t>(t_idx)][b_idx] = curr;
      curr += local_counts[static_cast<size_t>(t_idx)][b_idx];
    }
  }

#pragma omp parallel default(none) shared(src, dest, offsets, n, pass, kBits)
  {
    int tid = omp_get_thread_num();
#pragma omp for
    for (size_t i = 0; i < n; ++i) {
      size_t bucket = (src[i] >> (static_cast<size_t>(pass) * kBits)) & 255;
      dest[offsets[static_cast<size_t>(tid)][bucket]++] = src[i];
    }
  }
}
}  // namespace

TitaevSortirovkaBetcheraOMP::TitaevSortirovkaBetcheraOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool TitaevSortirovkaBetcheraOMP::ValidationImpl() {
  return !GetInput().empty();
}

bool TitaevSortirovkaBetcheraOMP::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

void TitaevSortirovkaBetcheraOMP::ConvertToKeys(const InType &input, std::vector<uint64_t> &keys) {
  const size_t n = input.size();
#pragma omp parallel for default(none) shared(keys, input, n)
  for (size_t i = 0; i < n; i++) {
    keys[i] = DoubleToOrderedUint(input[i]);
  }
}

void TitaevSortirovkaBetcheraOMP::RadixSortParallel(std::vector<uint64_t> &keys) {
  const size_t n = keys.size();
  if (n <= 1) {
    return;
  }
  std::vector<uint64_t> tmp(n);
  for (int pass = 0; pass < 8; ++pass) {
    if (pass % 2 == 0) {
      RadixPass(pass, n, keys, tmp);
    } else {
      RadixPass(pass, n, tmp, keys);
    }
  }
}

void TitaevSortirovkaBetcheraOMP::ConvertFromKeys(const std::vector<uint64_t> &keys, OutType &output) {
  const size_t n = keys.size();
  output.resize(n);
#pragma omp parallel for default(none) shared(keys, output, n)
  for (size_t i = 0; i < n; ++i) {
    output[i] = OrderedUintToDouble(keys[i]);
  }
}

void TitaevSortirovkaBetcheraOMP::BatcherStepParallel(OutType &result, size_t n, size_t step, size_t stage) {
#pragma omp parallel for default(none) shared(result, n, step, stage)
  for (size_t i = 0; i < n; ++i) {
    size_t j = i ^ stage;
    if (j > i && j < n) {
      const bool ascending = (i & step) == 0;
      if (ascending ? (result[i] > result[j]) : (result[i] < result[j])) {
        std::swap(result[i], result[j]);
      }
    }
  }
}

void TitaevSortirovkaBetcheraOMP::BatcherSortParallel() {
  auto &result = GetOutput();
  const size_t n = result.size();
  for (size_t step = 1; step < n; step <<= 1) {
    for (size_t stage = step; stage > 0; stage >>= 1) {
      BatcherStepParallel(result, n, step, stage);
    }
  }
}

bool TitaevSortirovkaBetcheraOMP::RunImpl() {
  auto &input = GetInput();
  const size_t n = input.size();
  if (n <= 1) {
    return true;
  }

  std::vector<uint64_t> keys(n);
  ConvertToKeys(input, keys);
  RadixSortParallel(keys);
  ConvertFromKeys(keys, GetOutput());

  if ((n & (n - 1)) == 0) {
    BatcherSortParallel();
  }
  return true;
}

bool TitaevSortirovkaBetcheraOMP::PostProcessingImpl() {
  return true;
}

}  // namespace titaev_m_sortirovka_betchera
