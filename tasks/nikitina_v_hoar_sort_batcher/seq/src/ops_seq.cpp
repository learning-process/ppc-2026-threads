#include "nikitina_v_hoar_sort_batcher/seq/include/ops_seq.hpp"

#include <algorithm>
#include <utility>
#include <vector>

#include "nikitina_v_hoar_sort_batcher/common/include/common.hpp"

namespace nikitina_v_hoar_sort_batcher {

namespace {

void QuickSortHoare(std::vector<int> &arr, int low, int high) {
  std::vector<std::pair<int, int>> stack;
  stack.emplace_back(low, high);

  while (!stack.empty()) {
    auto [l, h] = stack.back();
    stack.pop_back();

    if (l >= h) {
      continue;
    }

    int pivot = arr[l + ((h - l) / 2)];
    int i = l - 1;
    int j = h + 1;

    while (true) {
      i++;
      while (arr[i] < pivot) {
        i++;
      }

      j--;
      while (arr[j] > pivot) {
        j--;
      }

      if (i >= j) {
        break;
      }
      std::swap(arr[i], arr[j]);
    }

    stack.emplace_back(l, j);
    stack.emplace_back(j + 1, h);
  }
}

}  // namespace

HoareSortBatcherSEQ::HoareSortBatcherSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool HoareSortBatcherSEQ::ValidationImpl() {
  return true;
}

bool HoareSortBatcherSEQ::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool HoareSortBatcherSEQ::RunImpl() {
  auto &out = GetOutput();
  if (!out.empty()) {
    QuickSortHoare(out, 0, static_cast<int>(out.size()) - 1);
  }
  return true;
}

bool HoareSortBatcherSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace nikitina_v_hoar_sort_batcher
