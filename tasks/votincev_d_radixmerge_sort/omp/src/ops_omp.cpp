#include "votincev_d_radixmerge_sort/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "votincev_d_radixmerge_sort/common/include/common.hpp"

namespace votincev_d_radixmerge_sort {

VotincevDRadixMergeSortOMP::VotincevDRadixMergeSortOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool VotincevDRadixMergeSortOMP::ValidationImpl() {
  return !GetInput().empty();
}

bool VotincevDRadixMergeSortOMP::PreProcessingImpl() {
  return true;
}

void VotincevDRadixMergeSortOMP::SortByDigit(std::vector<int32_t> &array, int32_t exp) {
  std::vector<std::vector<int32_t>> buckets(10);

  for (const auto &num : array) {
    int32_t digit = (num / exp) % 10;
    buckets[digit].push_back(num);
  }

  size_t index = 0;
  for (int i = 0; i < 10; ++i) {
    for (const auto &val : buckets[i]) {
      array[index++] = val;
    }
    buckets[i].clear();
  }
}

bool VotincevDRadixMergeSortOMP::RunImpl() {
  if (GetInput().empty()) {
    return false;
  }

  std::vector<int32_t> working_array = GetInput();
  auto n = static_cast<int32_t>(working_array.size());

  int32_t min_val = working_array[0];
#pragma omp parallel for reduction(min : min_val) default(none) shared(working_array, n)
  for (int32_t i = 0; i < n; ++i) {
    min_val = std::min(min_val, working_array[i]);
  }

  if (min_val < 0) {
#pragma omp parallel for default(none) shared(working_array, n, min_val)
    for (int32_t i = 0; i < n; ++i) {
      working_array[i] -= min_val;
    }
  }

  int32_t max_val = working_array[0];
#pragma omp parallel for reduction(max : max_val) default(none) shared(working_array, n)
  for (int32_t i = 0; i < n; ++i) {
    max_val = std::max(max_val, working_array[i]);
  }

  for (int32_t exp = 1; max_val / exp > 0; exp *= 10) {
    SortByDigit(working_array, exp);
  }

  if (min_val < 0) {
#pragma omp parallel for default(none) shared(working_array, n, min_val)
    for (int32_t i = 0; i < n; ++i) {
      working_array[i] += min_val;
    }
  }

  GetOutput() = working_array;

  return true;
}

bool VotincevDRadixMergeSortOMP::PostProcessingImpl() {
  return true;
}

}  // namespace votincev_d_radixmerge_sort
