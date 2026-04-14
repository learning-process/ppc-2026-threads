#include "frolova_s_radix_sort_double/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "frolova_s_radix_sort_double/common/include/common.hpp"

namespace frolova_s_radix_sort_double {

namespace {

void RadixSortChunk(std::vector<double> &chunk) {
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

void ProcessChunk(std::vector<double> &chunk) {
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

}  // namespace

FrolovaSRadixSortDoubleOMP::FrolovaSRadixSortDoubleOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool FrolovaSRadixSortDoubleOMP::ValidationImpl() {
  return !GetInput().empty();
}

bool FrolovaSRadixSortDoubleOMP::PreProcessingImpl() {
  return true;
}

bool FrolovaSRadixSortDoubleOMP::RunImpl() {
  const std::vector<double> &input = GetInput();
  if (input.empty()) {
    return false;
  }

  size_t n = input.size();
  std::vector<double> working = input;

  int max_threads = omp_get_max_threads();
  int num_threads_to_use = std::min(max_threads, std::max(1, static_cast<int>(n / 10000)));
  if (num_threads_to_use == 0) {
    num_threads_to_use = 1;
  }

  std::vector<size_t> chunk_sizes(num_threads_to_use, n / num_threads_to_use);
  for (size_t i = 0; i < n % static_cast<size_t>(num_threads_to_use); i++) {
    chunk_sizes[i]++;
  }

  std::vector<size_t> chunk_offsets(num_threads_to_use, 0);
  for (int i = 1; i < num_threads_to_use; i++) {
    chunk_offsets[i] = chunk_offsets[i - 1] + chunk_sizes[i - 1];
  }

#pragma omp parallel num_threads(num_threads_to_use) default(none) \
    shared(working, chunk_sizes, chunk_offsets, num_threads_to_use)
  {
    int tid = omp_get_thread_num();
    if (tid < num_threads_to_use) {
      size_t offset = chunk_offsets[tid];
      size_t size = chunk_sizes[tid];

      std::vector<double> chunk(working.begin() + static_cast<ptrdiff_t>(offset),
                                working.begin() + static_cast<ptrdiff_t>(offset + size));

      ProcessChunk(chunk);

      for (size_t i = 0; i < size; ++i) {
        working[offset + i] = chunk[i];
      }
    }
  }

  std::vector<double> merged_result;
  merged_result.assign(working.begin(), working.begin() + static_cast<ptrdiff_t>(chunk_sizes[0]));

  for (int i = 1; i < num_threads_to_use; i++) {
    std::vector<double> merged(merged_result.size() + chunk_sizes[i]);
    auto next_chunk_begin = working.begin() + static_cast<ptrdiff_t>(chunk_offsets[i]);
    auto next_chunk_end = next_chunk_begin + static_cast<ptrdiff_t>(chunk_sizes[i]);
    std::merge(merged_result.begin(), merged_result.end(), next_chunk_begin, next_chunk_end, merged.begin());

    merged_result = std::move(merged);
  }

  GetOutput() = std::move(merged_result);
  return true;
}

bool FrolovaSRadixSortDoubleOMP::PostProcessingImpl() {
  return true;
}

}  // namespace frolova_s_radix_sort_double
