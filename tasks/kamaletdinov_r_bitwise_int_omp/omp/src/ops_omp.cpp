#include "kamaletdinov_r_bitwise_int_omp/omp/include/ops_omp.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <utility>
#include <vector>

#include "kamaletdinov_r_bitwise_int_omp/common/include/common.hpp"
#include "util/include/util.hpp"

namespace kamaletdinov_r_bitwise_int_omp {

namespace {

void CountingSortByDigit(std::vector<int> &data, int exp) {
  int size = static_cast<int>(data.size());
  std::vector<int> output(size);
  std::array<int, 10> count = {};

  for (int idx = 0; idx < size; idx++) {
    count.at((data[idx] / exp) % 10)++;
  }

  for (int idx = 1; idx < 10; idx++) {
    count.at(idx) += count.at(idx - 1);
  }

  for (int idx = size - 1; idx >= 0; idx--) {
    int digit = (data[idx] / exp) % 10;
    output[count.at(digit) - 1] = data[idx];
    count.at(digit)--;
  }

  data = output;
}

void RadixSortPositive(std::vector<int> &data) {
  if (data.empty()) {
    return;
  }

  int max_value = *std::ranges::max_element(data);

  for (int exp = 1; max_value / exp > 0; exp *= 10) {
    CountingSortByDigit(data, exp);
  }
}

void RadixSortChunk(std::vector<int> &chunk) {
  if (chunk.size() <= 1) {
    return;
  }

  std::vector<int> negative;
  std::vector<int> positive;

  negative.reserve(chunk.size());
  positive.reserve(chunk.size());

  for (int value : chunk) {
    if (value < 0) {
      negative.push_back(-value);
    } else {
      positive.push_back(value);
    }
  }

  RadixSortPositive(negative);
  RadixSortPositive(positive);

  std::size_t index = 0;
  for (int idx = static_cast<int>(negative.size()) - 1; idx >= 0; idx--) {
    chunk[index++] = -negative[static_cast<std::size_t>(idx)];
  }
  for (int value : positive) {
    chunk[index++] = value;
  }
}

std::vector<int> MergeSorted(const std::vector<int> &left, const std::vector<int> &right) {
  std::vector<int> result;
  result.reserve(left.size() + right.size());
  std::size_t li = 0;
  std::size_t ri = 0;
  while (li < left.size() && ri < right.size()) {
    if (left[li] <= right[ri]) {
      result.push_back(left[li++]);
    } else {
      result.push_back(right[ri++]);
    }
  }
  while (li < left.size()) {
    result.push_back(left[li++]);
  }
  while (ri < right.size()) {
    result.push_back(right[ri++]);
  }
  return result;
}

}  // namespace

void BitwiseSort(std::vector<int> &data) {
  if (data.size() <= 1) {
    return;
  }

  int num_threads = ppc::util::GetNumThreads();
  int total = static_cast<int>(data.size());

  if (num_threads <= 1 || total < num_threads) {
    RadixSortChunk(data);
    return;
  }

  std::vector<std::vector<int>> chunks(num_threads);
  int base_size = total / num_threads;
  int remainder = total % num_threads;

  int offset = 0;
  for (int ti = 0; ti < num_threads; ti++) {
    int chunk_size = base_size + (ti < remainder ? 1 : 0);
    chunks[ti].assign(data.begin() + offset, data.begin() + offset + chunk_size);
    offset += chunk_size;
  }

#pragma omp parallel for schedule(static) default(none) shared(chunks, num_threads)
  for (int ti = 0; ti < num_threads; ti++) {
    RadixSortChunk(chunks[ti]);
  }

  std::vector<int> merged = std::move(chunks[0]);
  for (int ti = 1; ti < num_threads; ti++) {
    merged = MergeSorted(merged, chunks[ti]);
  }

  data = std::move(merged);
}

KamaletdinovRBitwiseIntOMP::KamaletdinovRBitwiseIntOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool KamaletdinovRBitwiseIntOMP::ValidationImpl() {
  return GetInput() >= 0;
}

bool KamaletdinovRBitwiseIntOMP::PreProcessingImpl() {
  int size = GetInput();
  data_.clear();
  if (size <= 0) {
    return true;
  }
  data_.reserve(static_cast<std::size_t>(size));
  for (int idx = 0; idx < size; idx++) {
    data_.push_back(size - idx);
  }
  return true;
}

bool KamaletdinovRBitwiseIntOMP::RunImpl() {
  if (!data_.empty()) {
    BitwiseSort(data_);
  }
  return true;
}

bool KamaletdinovRBitwiseIntOMP::PostProcessingImpl() {
  bool sorted = std::ranges::is_sorted(data_);
  GetOutput() = sorted ? GetInput() : 0;
  return true;
}

}  // namespace kamaletdinov_r_bitwise_int_omp
