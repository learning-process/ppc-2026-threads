#include "spichek_d_radix_sort_for_integers_with_simple_merging/seq/include/ops_seq.hpp"

#include <algorithm>
#include <vector>

#include "spichek_d_radix_sort_for_integers_with_simple_merging/common/include/common.hpp"

namespace spichek_d_radix_sort_for_integers_with_simple_merging {

SpichekDRadixSortSEQ::SpichekDRadixSortSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool SpichekDRadixSortSEQ::ValidationImpl() {
  return !GetInput().empty();
}

bool SpichekDRadixSortSEQ::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool SpichekDRadixSortSEQ::RunImpl() {
  if (GetOutput().empty()) {
    return true;
  }

  RadixSort(GetOutput());
  return true;
}

bool SpichekDRadixSortSEQ::PostProcessingImpl() {
  return std::ranges::is_sorted(GetOutput());
}

void SpichekDRadixSortSEQ::RadixSort(std::vector<int> &data) {
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

  for (int exp = 1; max_val / exp > 0; exp *= 10) {
    std::vector<int> output(data.size());
    std::vector<int> count(10, 0);

    for (int x : data) {
      count[(x / exp) % 10]++;
    }

    for (int i = 1; i < 10; i++) {
      count[i] += count[i - 1];
    }

    for (int i = static_cast<int>(data.size()) - 1; i >= 0; i--) {
      output[count[(data[i] / exp) % 10] - 1] = data[i];
      count[(data[i] / exp) % 10]--;
    }

    data = output;
  }

  if (min_val < 0) {
    for (auto &x : data) {
      x += min_val;
    }
  }
}

}  // namespace spichek_d_radix_sort_for_integers_with_simple_merging
