#include "kamaletdinov_r_bitwise_int_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>

#include "kamaletdinov_r_bitwise_int_seq/common/include/common.hpp"

namespace kamaletdinov_r_bitwise_int_seq {

namespace {

void CountingSortByDigit(std::vector<int> &arr, int exp) {
  int n = static_cast<int>(arr.size());
  std::vector<int> output(n);
  std::array<int, 10> count = {};

  for (int i = 0; i < n; i++) {
    count.at((arr[i] / exp) % 10)++;
  }

  for (int i = 1; i < 10; i++) {
    count.at(i) += count.at(i - 1);
  }

  for (int i = n - 1; i >= 0; i--) {
    output[count.at((arr[i] / exp) % 10) - 1] = arr[i];
    count.at((arr[i] / exp) % 10)--;
  }

  arr = output;
}

void RadixSortPositive(std::vector<int> &arr) {
  if (arr.empty()) {
    return;
  }

  int max_val = *std::ranges::max_element(arr);

  for (int exp = 1; max_val / exp > 0; exp *= 10) {
    CountingSortByDigit(arr, exp);
    if (exp > max_val / 10) {
      break;
    }
  }
}

}  // namespace

void BitwiseSort(std::vector<int> &arr) {
  if (arr.size() <= 1) {
    return;
  }

  std::vector<int> neg;
  std::vector<int> pos;

  for (int x : arr) {
    if (x < 0) {
      neg.push_back(-x);
    } else {
      pos.push_back(x);
    }
  }

  RadixSortPositive(neg);
  RadixSortPositive(pos);

  std::size_t idx = 0;
  for (int i = static_cast<int>(neg.size()) - 1; i >= 0; i--) {
    arr[idx++] = -neg[i];
  }
  for (int p : pos) {
    arr[idx++] = p;
  }
}

KamaletdinovRBitwiseIntSEQ::KamaletdinovRBitwiseIntSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool KamaletdinovRBitwiseIntSEQ::ValidationImpl() {
  return GetInput() >= 0;
}

bool KamaletdinovRBitwiseIntSEQ::PreProcessingImpl() {
  int n = GetInput();
  data_.resize(n);
  for (int i = 0; i < n; i++) {
    data_[i] = (n / 2) - i;
  }
  return true;
}

bool KamaletdinovRBitwiseIntSEQ::RunImpl() {
  BitwiseSort(data_);
  return true;
}

bool KamaletdinovRBitwiseIntSEQ::PostProcessingImpl() {
  bool sorted = std::ranges::is_sorted(data_);
  GetOutput() = sorted ? GetInput() : 0;
  return true;
}

}  // namespace kamaletdinov_r_bitwise_int_seq
