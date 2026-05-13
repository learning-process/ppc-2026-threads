#include "kamaletdinov_r_bitwise_int/tbb/include/ops_tbb.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>

#include "kamaletdinov_r_bitwise_int/common/include/common.hpp"
#include "oneapi/tbb/parallel_for.h"
#include "oneapi/tbb/task_arena.h"

namespace kamaletdinov_r_bitwise_int {

namespace {

void CountingSortByDigit(std::vector<int> &arr, int exp) {
  const int n = static_cast<int>(arr.size());
  const int num_parts = tbb::this_task_arena::max_concurrency();

  std::vector<std::array<int, 10>> part_counts(num_parts);

  tbb::parallel_for(0, num_parts, [&](int part) {
    part_counts.at(part).fill(0);
    const int begin = part * n / num_parts;
    const int end = (part + 1) * n / num_parts;
    for (int i = begin; i < end; i++) {
      const int digit = (arr.at(i) / exp) % 10;
      part_counts.at(part).at(digit)++;
    }
  });

  std::array<int, 10> total = {};
  for (int pi = 0; pi < num_parts; pi++) {
    for (int di = 0; di < 10; di++) {
      total.at(di) += part_counts.at(pi).at(di);
    }
  }

  std::array<int, 10> starts = {};
  for (int di = 1; di < 10; di++) {
    starts.at(di) = starts.at(di - 1) + total.at(di - 1);
  }

  std::vector<std::array<int, 10>> part_positions(num_parts);
  for (int di = 0; di < 10; di++) {
    int off = starts.at(di);
    for (int pi = 0; pi < num_parts; pi++) {
      part_positions.at(pi).at(di) = off;
      off += part_counts.at(pi).at(di);
    }
  }

  std::vector<int> output(n);

  tbb::parallel_for(0, num_parts, [&](int part) {
    auto pos = part_positions.at(part);
    const int begin = part * n / num_parts;
    const int end = (part + 1) * n / num_parts;
    for (int i = begin; i < end; i++) {
      const int digit = (arr.at(i) / exp) % 10;
      output.at(pos.at(digit)++) = arr.at(i);
    }
  });

  arr.swap(output);
}

void RadixSortPositive(std::vector<int> &arr) {
  if (arr.empty()) {
    return;
  }

  const int max_val = *std::ranges::max_element(arr);
  for (int exp = 1; max_val / exp > 0; exp *= 10) {
    CountingSortByDigit(arr, exp);
    if (exp > max_val / 10) {
      break;
    }
  }
}

void BitwiseSortTBB(std::vector<int> &arr) {
  if (arr.size() <= 1) {
    return;
  }

  std::vector<int> neg;
  std::vector<int> pos;
  neg.reserve(arr.size());
  pos.reserve(arr.size());

  for (int val : arr) {
    if (val < 0) {
      neg.push_back(-val);
    } else {
      pos.push_back(val);
    }
  }

  RadixSortPositive(neg);
  RadixSortPositive(pos);

  std::size_t idx = 0;
  for (int i = static_cast<int>(neg.size()) - 1; i >= 0; i--) {
    arr.at(idx++) = -neg.at(i);
  }
  for (int v : pos) {
    arr.at(idx++) = v;
  }
}

}  // namespace

KamaletdinovRBitwiseIntTBB::KamaletdinovRBitwiseIntTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool KamaletdinovRBitwiseIntTBB::ValidationImpl() {
  return GetInput() >= 0;
}

bool KamaletdinovRBitwiseIntTBB::PreProcessingImpl() {
  const int n = GetInput();
  data_.resize(n);
  tbb::parallel_for(0, n, [&](int i) { data_.at(i) = (n / 2) - i; });
  return true;
}

bool KamaletdinovRBitwiseIntTBB::RunImpl() {
  BitwiseSortTBB(data_);
  return true;
}

bool KamaletdinovRBitwiseIntTBB::PostProcessingImpl() {
  const bool sorted = std::ranges::is_sorted(data_);
  GetOutput() = sorted ? GetInput() : 0;
  return true;
}

}  // namespace kamaletdinov_r_bitwise_int
