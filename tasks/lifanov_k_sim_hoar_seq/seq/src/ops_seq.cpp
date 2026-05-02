#include "lifanov_k_simple_hoar_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <vector>

#include "lifanov_k_simple_hoar_seq/common/include/common.hpp"

namespace lifanov_k_simple_hoar_seq {

LifanovKSimpleHoarSEQ::LifanovKSimpleHoarSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool LifanovKSimpleHoarSEQ::ValidationImpl() {
  return !GetInput().empty();
}

bool LifanovKSimpleHoarSEQ::PreProcessingImpl() {
  data_ = GetInput();
  return true;
}

int LifanovKSimpleHoarSEQ::Partition(std::vector<int> &arr, int low, int high) {
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
      return j;
    }
    std::swap(arr[i], arr[j]);
  }
}

void LifanovKSimpleHoarSEQ::QuickSortHoare(std::vector<int> &arr, int low, int high) {
  if (low < high) {
    int p = Partition(arr, low, high);
    QuickSortHoare(arr, low, p);
    QuickSortHoare(arr, p + 1, high);
  }
}

std::vector<int> LifanovKSimpleHoarSEQ::Merge(const std::vector<int> &left, const std::vector<int> &right) {
  std::vector<int> res;
  res.reserve(left.size() + right.size());
  size_t i = 0, j = 0;
  while (i < left.size() && j < right.size()) {
    if (left[i] <= right[j]) {
      res.push_back(left[i++]);
    } else {
      res.push_back(right[j++]);
    }
  }
  while (i < left.size()) {
    res.push_back(left[i++]);
  }
  while (j < right.size()) {
    res.push_back(right[j++]);
  }
  return res;
}

bool LifanovKSimpleHoarSEQ::RunImpl() {
  if (data_.size() <= 1) {
    return true;
  }

  size_t mid = data_.size() / 2;
  std::vector<int> left_part(data_.begin(), data_.begin() + mid);
  std::vector<int> right_part(data_.begin() + mid, data_.end());

  QuickSortHoare(left_part, 0, static_cast<int>(left_part.size() - 1));
  QuickSortHoare(right_part, 0, static_cast<int>(right_part.size() - 1));

  data_ = Merge(left_part, right_part);

  return true;
}

bool LifanovKSimpleHoarSEQ::PostProcessingImpl() {
  GetOutput() = data_;
  return std::is_sorted(GetOutput().begin(), GetOutput().end());
}

}  // namespace lifanov_k_simple_hoar_seq
