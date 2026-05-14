#include "titaev_m_sortirovka_betchera/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cstring>

namespace titaev_m_sortirovka_betchera {

namespace {
uint64_t DoubleToOrderedUint(double value) {
  uint64_t x;
  std::memcpy(&x, &value, sizeof(double));
  constexpr uint64_t kSignMask = (1ULL << 63);
  return ((x & kSignMask) != 0ULL) ? ~x : (x ^ kSignMask);
}

double OrderedUintToDouble(uint64_t x) {
  constexpr uint64_t kSignMask = (1ULL << 63);
  x = ((x & kSignMask) != 0ULL) ? (x ^ kSignMask) : ~x;
  double result;
  std::memcpy(&result, &x, sizeof(double));
  return result;
}
}  // namespace

bool TitaevSortirovkaBetcheraOMP::ValidationImpl() {
  return !task_data->inputs.empty() && !task_data->outputs.empty() &&
         task_data->inputs_count[0] == task_data->outputs_count[0];
}

bool TitaevSortirovkaBetcheraOMP::PreProcessingImpl() {
  auto *ptr = reinterpret_cast<double *>(task_data->inputs[0]);
  size_t n = task_data->inputs_count[0];
  GetInput().assign(ptr, ptr + n);
  GetOutput().resize(n);
  return true;
}

void TitaevSortirovkaBetcheraOMP::ConvertToKeys(const InType &input, std::vector<uint64_t> &keys) {
  size_t n = input.size();
#pragma omp parallel for
  for (int i = 0; i < (int)n; i++) {
    keys[i] = DoubleToOrderedUint(input[i]);
  }
}

void TitaevSortirovkaBetcheraOMP::RadixSort(std::vector<uint64_t> &keys) {
  size_t n = keys.size();
  std::vector<uint64_t> tmp(n);
  for (int pass = 0; pass < 8; pass++) {
    size_t count[256] = {0};
    int shift = pass * 8;
    for (size_t i = 0; i < n; i++) {
      count[(keys[i] >> shift) & 0xFF]++;
    }
    for (int i = 1; i < 256; i++) {
      count[i] += count[i - 1];
    }
    for (int i = (int)n - 1; i >= 0; i--) {
      tmp[--count[(keys[i] >> shift) & 0xFF]] = keys[i];
    }
    keys = tmp;
  }
}

void TitaevSortirovkaBetcheraOMP::ConvertFromKeys(const std::vector<uint64_t> &keys, OutType &output) {
#pragma omp parallel for
  for (int i = 0; i < (int)keys.size(); i++) {
    output[i] = OrderedUintToDouble(keys[i]);
  }
}

void TitaevSortirovkaBetcheraOMP::BatcherStep(OutType &result, size_t n, size_t step, size_t stage) {
#pragma omp parallel for
  for (int i = 0; i < (int)n; i++) {
    size_t j = i ^ stage;
    if (j > (size_t)i && j < n) {
      bool ascending = (i & step) == 0;
      if (ascending ? (result[i] > result[j]) : (result[i] < result[j])) {
        std::swap(result[i], result[j]);
      }
    }
  }
}

void TitaevSortirovkaBetcheraOMP::BatcherSort() {
  size_t n = GetOutput().size();
  for (size_t step = 1; step < n; step <<= 1) {
    for (size_t stage = step; stage > 0; stage >>= 1) {
      BatcherStep(GetOutput(), n, step, stage);
    }
  }
}

bool TitaevSortirovkaBetcheraOMP::RunImpl() {
  size_t n = GetInput().size();
  if (n == 0) {
    return true;
  }
  std::vector<uint64_t> keys(n);
  ConvertToKeys(GetInput(), keys);
  RadixSort(keys);
  ConvertFromKeys(keys, GetOutput());
  if ((n & (n - 1)) == 0) {
    BatcherSort();
  }
  return true;
}

bool TitaevSortirovkaBetcheraOMP::PostProcessingImpl() {
  auto *ptr = reinterpret_cast<double *>(task_data->outputs[0]);
  std::copy(GetOutput().begin(), GetOutput().end(), ptr);
  return true;
}

}  // namespace titaev_m_sortirovka_betchera
