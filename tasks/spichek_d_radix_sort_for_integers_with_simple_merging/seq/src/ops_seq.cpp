#include "spichek_d_radix_sort_for_integers_with_simple_merging/seq/include/ops_seq.hpp"

#include <algorithm>
#include <vector>

namespace spichek_d_radix_sort_for_integers_with_simple_merging {

RadixSortSimpleMergingSEQ::RadixSortSimpleMergingSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool RadixSortSimpleMergingSEQ::ValidationImpl() {
  return true;
}

bool RadixSortSimpleMergingSEQ::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool RadixSortSimpleMergingSEQ::RunImpl() {
  if (GetOutput().empty()) {
    return true;
  }
  auto &data = GetOutput();
  int min_val = *std::min_element(data.begin(), data.end());
  if (min_val < 0) {
    for (auto &x : data) {
      x -= min_val;
    }
  }
  int max_val = *std::max_element(data.begin(), data.end());

  for (int exp = 1; max_val / exp > 0; exp *= 10) {
    std::vector<int> output(data.size());
    int count[10] = {0};
    for (int i = 0; i < (int)data.size(); i++) {
      count[(data[i] / exp) % 10]++;
    }
    for (int i = 1; i < 10; i++) {
      count[i] += count[i - 1];
    }
    for (int i = (int)data.size() - 1; i >= 0; i--) {
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
  return true;
}

bool RadixSortSimpleMergingSEQ::PostProcessingImpl() {
  return std::is_sorted(GetOutput().begin(), GetOutput().end());
}

}  // namespace spichek_d_radix_sort_for_integers_with_simple_merging
