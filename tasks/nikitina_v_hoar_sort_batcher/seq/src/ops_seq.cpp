#include "nikitina_v_hoar_sort_batcher/seq/include/ops_seq.hpp"

#include <algorithm>
#include <utility>

#include "nikitina_v_hoar_sort_batcher/common/include/common.hpp"

namespace nikitina_v_hoar_sort_batcher {

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

void HoareSortBatcherSEQ::QuickSortHoare(std::vector<int> &arr, int low, int high) {
  if (low >= high) {
    return;
  }

  int pivot = arr[low + (high - low) / 2];
  int i = low - 1;
  int j = high + 1;

  while (true) {
    do {
      i++;
    } while (arr[i] < pivot);
    do {
      j--;
    } while (arr[j] > pivot);

    if (i >= j) {
      break;
    }
    std::swap(arr[i], arr[j]);
  }

  QuickSortHoare(arr, low, j);
  QuickSortHoare(arr, j + 1, high);
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
