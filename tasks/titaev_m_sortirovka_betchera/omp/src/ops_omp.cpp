#include "titaev_m_sortirovka_betchera/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cstring>
#include <limits>
#include <vector>

namespace titaev_m_sortirovka_betchera {

TitaevSortirovkaBetcheraOMP::TitaevSortirovkaBetcheraOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
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
  uint64_t bits;
  std::memcpy(&bits, &v, 8);
  return (bits & 0x8000000000000000ULL) ? ~bits : (bits | 0x8000000000000000ULL);
}

double TitaevSortirovkaBetcheraOMP::UnpackDouble(uint64_t k) noexcept {
  uint64_t bits = (k & 0x8000000000000000ULL) ? (k & ~0x8000000000000000ULL) : ~k;
  double v;
  std::memcpy(&v, &bits, 8);
  return v;
}

void TitaevSortirovkaBetcheraOMP::CompareSwap(std::vector<double> &arr, size_t i, size_t j) {
  if (i < arr.size() && j < arr.size() && arr[i] > arr[j]) {
    std::swap(arr[i], arr[j]);
  }
}

void TitaevSortirovkaBetcheraOMP::LSDRadixSort(std::vector<double> &array) {
  const size_t n = array.size();
  if (n <= 1) {
    return;
  }

  std::vector<uint64_t> keys(n);
  for (size_t i = 0; i < n; ++i) {
    keys[i] = PackDouble(array[i]);
  }

  std::vector<uint64_t> tmp_keys(n);
  std::vector<double> tmp_vals(n);

  for (int pass = 0; pass < 8; ++pass) {
    size_t cnt[257] = {0};
    int shift = pass * 8;
    for (size_t i = 0; i < n; ++i) {
      cnt[((keys[i] >> shift) & 0xFF) + 1]++;
    }
    for (int i = 0; i < 256; ++i) {
      cnt[i + 1] += cnt[i];
    }
    for (size_t i = 0; i < n; ++i) {
      size_t pos = cnt[(keys[i] >> shift) & 0xFF]++;
      tmp_keys[pos] = keys[i];
      tmp_vals[pos] = array[i];
    }
    keys.swap(tmp_keys);
    array.swap(tmp_vals);
  }

  for (size_t i = 0; i < n; ++i) {
    array[i] = UnpackDouble(keys[i]);
  }
}

void TitaevSortirovkaBetcheraOMP::BatcherOddEvenMerge(std::vector<double> &arr, size_t n) {
  if (n <= 1) {
    return;
  }
  for (size_t po = n / 2; po > 0; po >>= 1) {
    if (po == n / 2) {
#pragma omp parallel for default(none) shared(arr, po)
      for (long long i = 0; i < (long long)po; ++i) {
        CompareSwap(arr, (size_t)i, (size_t)i + po);
      }
    } else {
#pragma omp parallel for default(none) shared(arr, n, po)
      for (long long i = (long long)po; i < (long long)(n - po); i += (long long)(2 * po)) {
        for (size_t j = 0; j < po; ++j) {
          CompareSwap(arr, (size_t)i + j, (size_t)i + j + po);
        }
      }
    }
  }
}

bool TitaevSortirovkaBetcheraOMP::RunImpl() {
  const auto &input_data = GetInput();
  if (input_data.empty()) {
    GetOutput() = std::vector<double>();
    return true;
  }

  size_t original_size = input_data.size();
  size_t pow2 = 1;
  while (pow2 < original_size) {
    pow2 <<= 1;
  }

  std::vector<double> data = input_data;
  data.resize(pow2, std::numeric_limits<double>::max());

  size_t half = pow2 / 2;
  if (half > 0) {
#pragma omp parallel sections shared(data)
    {
#pragma omp section
      {
        std::vector<double> l(data.begin(), data.begin() + half);
        LSDRadixSort(l);
        if (!l.empty()) {
          std::copy(l.begin(), l.end(), data.begin());
        }
      }
#pragma omp section
      {
        std::vector<double> r(data.begin() + half, data.end());
        LSDRadixSort(r);
        if (!r.empty()) {
          std::copy(r.begin(), r.end(), data.begin() + half);
        }
      }
    }
    BatcherOddEvenMerge(data, pow2);
  } else {
    LSDRadixSort(data);
  }

  data.resize(original_size);
  GetOutput() = std::move(data);
  return true;
}

}  // namespace titaev_m_sortirovka_betchera
