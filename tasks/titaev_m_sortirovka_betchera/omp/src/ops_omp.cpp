#include "titaev_m_sortirovka_betchera/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>

namespace titaev_m_sortirovka_betchera {

TitaevSortirovkaBetcheraOMP::TitaevSortirovkaBetcheraOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool TitaevSortirovkaBetcheraOMP::ValidationImpl() {
  return true;
}

bool TitaevSortirovkaBetcheraOMP::PreProcessingImpl() {
  return true;
}

bool TitaevSortirovkaBetcheraOMP::PostProcessingImpl() {
  return true;
}

uint64_t TitaevSortirovkaBetcheraOMP::PackDouble(double v) noexcept {
  uint64_t bits = 0ULL;
  std::memcpy(&bits, &v, sizeof(bits));
  if (bits & (1ULL << 63)) {
    bits = ~bits;
  } else {
    bits ^= (1ULL << 63);
  }
  return bits;
}

double TitaevSortirovkaBetcheraOMP::UnpackDouble(uint64_t k) noexcept {
  if (k & (1ULL << 63)) {
    k ^= (1ULL << 63);
  } else {
    k = ~k;
  }
  double v = 0.0;
  std::memcpy(&v, &k, sizeof(v));
  return v;
}

void TitaevSortirovkaBetcheraOMP::LSDRadixSort(std::vector<double> &array) {
  const size_t n = array.size();
  if (n <= 1) {
    return;
  }

  constexpr int kBits = 8;
  constexpr int kBuckets = 1 << kBits;
  constexpr int kPasses = static_cast<int>((sizeof(uint64_t) * 8) / kBits);

  std::vector<uint64_t> keys(n);
  for (size_t i = 0; i < n; ++i) {
    keys[i] = TitaevSortirovkaBetcheraOMP::PackDouble(array[i]);
  }

  std::vector<uint64_t> tmp_keys(n);
  std::vector<double> tmp_vals(n);

  for (int pass = 0; pass < kPasses; ++pass) {
    const int shift = pass * kBits;
    size_t cnt[kBuckets + 1] = {0};

    for (size_t i = 0; i < n; ++i) {
      size_t d = (keys[i] >> shift) & (kBuckets - 1);
      ++cnt[d + 1];
    }
    for (int i = 0; i < kBuckets; ++i) {
      cnt[i + 1] += cnt[i];
    }

    for (size_t i = 0; i < n; ++i) {
      size_t d = (keys[i] >> shift) & (kBuckets - 1);
      size_t pos = cnt[d]++;
      tmp_keys[pos] = keys[i];
      tmp_vals[pos] = array[i];
    }
    keys.swap(tmp_keys);
    array.swap(tmp_vals);
  }

  for (size_t i = 0; i < n; ++i) {
    array[i] = TitaevSortirovkaBetcheraOMP::UnpackDouble(keys[i]);
  }
}

void TitaevSortirovkaBetcheraOMP::BatcherOddEvenMerge(std::vector<double> &arr, size_t n) {
  for (size_t po = n / 2; po > 0; po >>= 1) {
    if (po == n / 2) {
#pragma omp parallel for default(none) shared(arr, po)
      for (ptrdiff_t i = 0; i < (ptrdiff_t)po; ++i) {
        CompareSwap(arr, (size_t)i, (size_t)i + po);
      }
    } else {
#pragma omp parallel for default(none) shared(arr, n, po)
      for (ptrdiff_t i = (ptrdiff_t)po; i < (ptrdiff_t)(n - po); i += (ptrdiff_t)(2 * po)) {
        for (size_t j = 0; j < po; ++j) {
          CompareSwap(arr, (size_t)i + j, (size_t)i + j + po);
        }
      }
    }
  }
}

bool TitaevSortirovkaBetcheraOMP::RunImpl() {
  auto data = GetInput();
  const size_t original_size = data.size();

  if (original_size <= 1) {
    GetOutput() = data;
    return true;
  }

  size_t pow2 = 1;
  while (pow2 < original_size) {
    pow2 <<= 1;
  }
  data.resize(pow2, std::numeric_limits<double>::max());

  const size_t half = pow2 / 2;
#pragma omp parallel sections default(none) shared(data, half)
  {
#pragma omp section
    {
      std::vector<double> left(data.begin(), data.begin() + half);
      LSDRadixSort(left);
      std::copy(left.begin(), left.end(), data.begin());
    }
#pragma omp section
    {
      std::vector<double> right(data.begin() + half, data.end());
      LSDRadixSort(right);
      std::copy(right.begin(), right.end(), data.begin() + half);
    }
  }

  BatcherOddEvenMerge(data, pow2);

  data.resize(original_size);
  GetOutput() = data;
  return true;
}

}  // namespace titaev_m_sortirovka_betchera
