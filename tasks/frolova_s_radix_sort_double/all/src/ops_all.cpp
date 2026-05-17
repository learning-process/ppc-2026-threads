#include "frolova_s_radix_sort_double/all/include/ops_all.hpp"

#include <omp.h>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "frolova_s_radix_sort_double/common/include/common.hpp"
#include "task/include/task.hpp"

namespace frolova_s_radix_sort_double {

FrolovaSRadixSortDoubleALL::FrolovaSRadixSortDoubleALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool FrolovaSRadixSortDoubleALL::ValidationImpl() {
  return !GetInput().empty();
}

bool FrolovaSRadixSortDoubleALL::PreProcessingImpl() {
  return true;
}

void FrolovaSRadixSortDoubleALL::RadixSortChunk(std::vector<double> &chunk) {
  const int radix = 256;
  const int num_bits = 8;
  const int num_passes = sizeof(uint64_t);

  std::vector<double> temp(chunk.size());

  for (int pass = 0; pass < num_passes; pass++) {
    std::vector<int> count(radix, 0);
    for (double value : chunk) {
      auto bits = std::bit_cast<uint64_t>(value);
      int byte = static_cast<int>((bits >> (pass * num_bits)) & 0xFF);
      count[byte]++;
    }
    int total = 0;
    for (int i = 0; i < radix; i++) {
      int old = count[i];
      count[i] = total;
      total += old;
    }
    for (double value : chunk) {
      auto bits = std::bit_cast<uint64_t>(value);
      int byte = static_cast<int>((bits >> (pass * num_bits)) & 0xFF);
      temp[count[byte]++] = value;
    }
    chunk.swap(temp);
  }
}

void FrolovaSRadixSortDoubleALL::ProcessChunk(std::vector<double> &chunk) {
  RadixSortChunk(chunk);

  std::vector<double> negative;
  std::vector<double> positive;
  negative.reserve(chunk.size());
  positive.reserve(chunk.size());

  for (double val : chunk) {
    if ((std::bit_cast<uint64_t>(val) >> 63) != 0U) {
      negative.push_back(val);
    } else {
      positive.push_back(val);
    }
  }

  std::ranges::reverse(negative);

  size_t pos = 0;
  for (double val : negative) {
    chunk[pos++] = val;
  }
  for (double val : positive) {
    chunk[pos++] = val;
  }
}

std::vector<double> FrolovaSRadixSortDoubleALL::SimpleMerge(const std::vector<double> &a,
                                                            const std::vector<double> &b) {
  std::vector<double> res;
  res.reserve(a.size() + b.size());

  size_t i = 0;
  size_t j = 0;

  while (i < a.size() && j < b.size()) {
    if (a[i] <= b[j]) {
      res.push_back(a[i++]);
    } else {
      res.push_back(b[j++]);
    }
  }
  while (i < a.size()) {
    res.push_back(a[i++]);
  }
  while (j < b.size()) {
    res.push_back(b[j++]);
  }
  return res;
}

std::vector<double> FrolovaSRadixSortDoubleALL::ParallelMerge(std::vector<std::vector<double>> &chunks) {
  if (chunks.empty()) {
    return {};
  }

  while (chunks.size() > 1) {
    std::vector<std::vector<double>> next_chunks;
    next_chunks.resize((chunks.size() + 1) / 2);

    size_t half_size = chunks.size() / 2;

#pragma omp parallel for default(none) shared(chunks, next_chunks, half_size)
    for (size_t i = 0; i < half_size; ++i) {
      next_chunks[i] = SimpleMerge(chunks[(2 * i)], chunks[(2 * i) + 1]);
    }

    if (chunks.size() % 2 != 0) {
      next_chunks.back() = std::move(chunks.back());
    }

    chunks = std::move(next_chunks);
  }

  return std::move(chunks[0]);
}

bool FrolovaSRadixSortDoubleALL::RunImpl() {
  const std::vector<double> &input = GetInput();
  if (input.empty()) {
    return true;
  }

  int num_chunks = omp_get_max_threads();
  if (num_chunks <= 0) {
    num_chunks = 1;
  }

  int total_size = static_cast<int>(input.size());
  int chunk_size = total_size / num_chunks;
  int remainder = total_size % num_chunks;

  std::vector<std::vector<double>> chunks(num_chunks);
  int offset = 0;
  for (int i = 0; i < num_chunks; ++i) {
    int cur_size = chunk_size + (i < remainder ? 1 : 0);
    if (cur_size > 0) {
      chunks[i].assign(input.begin() + offset, input.begin() + offset + cur_size);
      offset += cur_size;
    }
  }

  std::erase_if(chunks, [](const std::vector<double> &c) { return c.empty(); });

  num_chunks = static_cast<int>(chunks.size());

  tbb::parallel_for(0, num_chunks, [&](int i) { ProcessChunk(chunks[i]); });

  std::vector<double> sorted = ParallelMerge(chunks);
  GetOutput() = std::move(sorted);

  return true;
}

bool FrolovaSRadixSortDoubleALL::PostProcessingImpl() {
  return true;
}

}  // namespace frolova_s_radix_sort_double
