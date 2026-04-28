#include "spichek_d_radix_sort_for_integers_with_simple_merging/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <iterator>
#include <utility>
#include <vector>

#include "spichek_d_radix_sort_for_integers_with_simple_merging/common/include/common.hpp"
#include "util/include/util.hpp"

namespace spichek_d_radix_sort_for_integers_with_simple_merging {

SpichekDRadixSortOMP::SpichekDRadixSortOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool SpichekDRadixSortOMP::ValidationImpl() {
  return !GetInput().empty();
}

bool SpichekDRadixSortOMP::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool SpichekDRadixSortOMP::RunImpl() {
  if (GetOutput().empty()) {
    return true;
  }

  int num_threads = ppc::util::GetNumThreads();
  int n = static_cast<int>(GetOutput().size());

  if (num_threads <= 1 || n < num_threads) {
    RadixSort(GetOutput());
    return true;
  }

  std::vector<std::vector<int>> local_data(num_threads);
  int base_size = n / num_threads;
  int remainder = n % num_threads;

  std::vector<int> counts(num_threads);
  std::vector<int> displs(num_threads, 0);

  for (int i = 0; i < num_threads; ++i) {
    counts[i] = base_size + (i < remainder ? 1 : 0);
    if (i > 0) {
      displs[i] = displs[i - 1] + counts[i - 1];
    }
    local_data[i] = std::vector<int>(GetOutput().begin() + displs[i], GetOutput().begin() + displs[i] + counts[i]);
  }

#pragma omp parallel for num_threads(num_threads) default(none) shared(local_data, num_threads)
  for (int i = 0; i < num_threads; ++i) {
    RadixSort(local_data[i]);
  }

  std::vector<int> result = std::move(local_data[0]);
  for (int i = 1; i < num_threads; ++i) {
    std::vector<int> temp;
    temp.reserve(result.size() + local_data[i].size());
    std::ranges::merge(result, local_data[i], std::back_inserter(temp));
    result = std::move(temp);
  }

  GetOutput() = std::move(result);
  return true;
}

bool SpichekDRadixSortOMP::PostProcessingImpl() {
  return std::ranges::is_sorted(GetOutput());
}

void SpichekDRadixSortOMP::RadixSort(std::vector<int> &data) {
  if (data.empty()) {
    return;
  }

  int min_val = *std::ranges::min_element(data);
  if (min_val < 0) {
    for (auto &x : data) {
      x -= min_val;
    }
  }

  int max_val = *std::ranges::max_element(data);

  for (int shift = 0; (max_val >> shift) > 0; shift += 8) {
    std::vector<int> output(data.size());
    std::vector<int> count(256, 0);

    for (int x : data) {
      count[(x >> shift) & 255]++;
    }

    for (int i = 1; i < 256; i++) {
      count[i] += count[i - 1];
    }

    for (int i = static_cast<int>(data.size()) - 1; i >= 0; i--) {
      output[count[(data[i] >> shift) & 255] - 1] = data[i];
      count[(data[i] >> shift) & 255]--;
    }

    data = std::move(output);
  }

  if (min_val < 0) {
    for (auto &x : data) {
      x += min_val;
    }
  }
}

}  // namespace spichek_d_radix_sort_for_integers_with_simple_merging
