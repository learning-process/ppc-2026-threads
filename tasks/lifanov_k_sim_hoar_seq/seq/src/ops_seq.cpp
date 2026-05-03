#include "lifanov_k_sim_hoar_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <vector>

#include "lifanov_k_sim_hoar_seq/common/include/common.hpp"

namespace lifanov_k_sim_hoar_seq {

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
  int mid = low + ((high - low) / 2);
  int pivot = arr[mid];

  int i = low - 1;
  int j = high + 1;

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
      return j;
    }

    std::swap(arr[i], arr[j]);
  }
}

void LifanovKSimpleHoarSEQ::QuickSortHoare(std::vector<int> &arr, int low, int high) {
  std::vector<std::pair<int, int>> stack;
  stack.push_back({low, high});

  while (!stack.empty()) {
    auto [l, r] = stack.back();
    stack.pop_back();

    if (l >= r) {
      continue;
    }

    int p = Partition(arr, l, r);

    if (l < p) {
      stack.push_back({l, p});
    }
    if (p + 1 < r) {
      stack.push_back({p + 1, r});
    }
  }
}

std::vector<int> LifanovKSimpleHoarSEQ::Merge(const std::vector<int> &left, const std::vector<int> &right) {
  std::vector<int> res;
  res.reserve(left.size() + right.size());

  size_t i = 0;
  size_t j = 0;

  while (i < left.size() && j < right.size()) {
    if (left[i] <= right[j]) {
      res.push_back(left[i]);
      i++;
    } else {
      res.push_back(right[j]);
      j++;
    }
  }

  while (i < left.size()) {
    res.push_back(left[i]);
    i++;
  }

  while (j < right.size()) {
    res.push_back(right[j]);
    j++;
  }

  return res;
}

bool LifanovKSimpleHoarSEQ::RunImpl() {
  if (data_.size() <= 1) {
    return true;
  }

  size_t mid = data_.size() / 2;
  auto mid_it = data_.begin() + static_cast<std::ptrdiff_t>(mid);

  std::vector<int> left_part(data_.begin(), mid_it);
  std::vector<int> right_part(mid_it, data_.end());

  QuickSortHoare(left_part, 0, static_cast<int>(left_part.size() - 1));
  QuickSortHoare(right_part, 0, static_cast<int>(right_part.size() - 1));

  data_ = Merge(left_part, right_part);

  return true;
}

bool LifanovKSimpleHoarSEQ::PostProcessingImpl() {
  GetOutput() = data_;
  return std::is_sorted(GetOutput().begin(), GetOutput().end());
}

}  // namespace lifanov_k_sim_hoar_seq
