#include "melnik_i_radix_sort_int/tbb/include/ops_tbb.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "melnik_i_radix_sort_int/common/include/common.hpp"
#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/parallel_for.h"
#include "util/include/util.hpp"

namespace melnik_i_radix_sort_int {

namespace {

constexpr int kBitsPerPass = 8;
constexpr std::size_t kBuckets = 1U << kBitsPerPass;

std::vector<MelnikIRadixSortIntTBB::Range> BuildInitialRanges(std::size_t data_size, int num_ranges) {
  std::vector<MelnikIRadixSortIntTBB::Range> ranges;
  ranges.reserve(static_cast<std::size_t>(num_ranges));

  for (int range_index = 0; range_index < num_ranges; ++range_index) {
    const std::size_t begin =
        (data_size * static_cast<std::size_t>(range_index)) / static_cast<std::size_t>(num_ranges);
    const std::size_t end =
        (data_size * (static_cast<std::size_t>(range_index) + 1U)) / static_cast<std::size_t>(num_ranges);
    ranges.push_back(MelnikIRadixSortIntTBB::Range{.begin = begin, .end = end});
  }

  return ranges;
}

}  // namespace

MelnikIRadixSortIntTBB::MelnikIRadixSortIntTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool MelnikIRadixSortIntTBB::ValidationImpl() {
  return !GetInput().empty();
}

bool MelnikIRadixSortIntTBB::PreProcessingImpl() {
  GetOutput() = GetInput();
  return !GetOutput().empty();
}

bool MelnikIRadixSortIntTBB::RunImpl() {
  if (GetOutput().empty()) {
    return false;
  }

  const std::size_t data_size = GetOutput().size();
  const int requested_threads = std::max(1, ppc::util::GetNumThreads());
  const int num_ranges = std::min(requested_threads, static_cast<int>(data_size));

  std::vector<int> buffer(data_size);
  auto &output = GetOutput();

  const std::vector<Range> ranges = BuildInitialRanges(data_size, num_ranges);

  tbb::parallel_for(tbb::blocked_range<std::size_t>(0, ranges.size()),
                    [&](const tbb::blocked_range<std::size_t> &range_block) {
    for (std::size_t range_index = range_block.begin(); range_index < range_block.end(); ++range_index) {
      const Range range = ranges[range_index];
      if (range.begin < range.end) {
        RadixSortRange(output, buffer, range.begin, range.end);
      }
    }
  });

  MergeSortedRanges(output, buffer, ranges);
  return !output.empty();
}

bool MelnikIRadixSortIntTBB::PostProcessingImpl() {
  return std::ranges::is_sorted(GetOutput());
}

void MelnikIRadixSortIntTBB::CountingSortByByte(const std::vector<int> &source, std::vector<int> &destination,
                                                std::size_t begin, std::size_t end, std::int64_t exp,
                                                std::int64_t offset) {
  std::array<std::size_t, kBuckets> count{};
  count.fill(0);

  for (std::size_t index = begin; index < end; ++index) {
    const std::int64_t shifted_value = static_cast<std::int64_t>(source[index]) + offset;
    const auto bucket = static_cast<std::size_t>((shifted_value / exp) % static_cast<std::int64_t>(kBuckets));
    ++count.at(bucket);
  }

  std::array<std::size_t, kBuckets> positions{};
  positions.at(0) = begin;
  for (std::size_t bucket = 1; bucket < kBuckets; ++bucket) {
    positions.at(bucket) = positions.at(bucket - 1U) + count.at(bucket - 1U);
  }

  for (std::size_t index = begin; index < end; ++index) {
    const std::int64_t shifted_value = static_cast<std::int64_t>(source[index]) + offset;
    const auto bucket = static_cast<std::size_t>((shifted_value / exp) % static_cast<std::int64_t>(kBuckets));
    destination[positions.at(bucket)] = source[index];
    ++positions.at(bucket);
  }
}

void MelnikIRadixSortIntTBB::RadixSortRange(std::vector<int> &data, std::vector<int> &buffer, std::size_t begin,
                                            std::size_t end) {
  if (end - begin <= 1) {
    return;
  }

  std::vector<int> *source = &data;
  std::vector<int> *destination = &buffer;
  const auto range_begin = data.begin() + static_cast<ptrdiff_t>(begin);
  const auto range_end = data.begin() + static_cast<ptrdiff_t>(end);
  const auto [min_it, max_it] = std::ranges::minmax_element(range_begin, range_end);
  const auto min_value = static_cast<std::int64_t>(*min_it);
  const auto max_value = static_cast<std::int64_t>(*max_it);
  const std::int64_t offset = (min_value < 0) ? -min_value : 0;
  const std::int64_t max_shifted_value = max_value + offset;

  for (std::int64_t exp = 1; max_shifted_value / exp > 0; exp <<= kBitsPerPass) {
    CountingSortByByte(*source, *destination, begin, end, exp, offset);
    std::swap(source, destination);
  }

  if (source != &data) {
    std::copy(source->begin() + static_cast<ptrdiff_t>(begin), source->begin() + static_cast<ptrdiff_t>(end),
              data.begin() + static_cast<ptrdiff_t>(begin));
  }
}

void MelnikIRadixSortIntTBB::MergeRanges(const std::vector<int> &source, std::vector<int> &destination, Range left,
                                         Range right, std::size_t write_begin) {
  std::size_t left_index = left.begin;
  std::size_t right_index = right.begin;
  std::size_t write_index = write_begin;

  while (left_index < left.end && right_index < right.end) {
    if (source[left_index] <= source[right_index]) {
      destination[write_index] = source[left_index];
      ++left_index;
    } else {
      destination[write_index] = source[right_index];
      ++right_index;
    }
    ++write_index;
  }

  if (left_index < left.end) {
    std::copy(source.begin() + static_cast<ptrdiff_t>(left_index), source.begin() + static_cast<ptrdiff_t>(left.end),
              destination.begin() + static_cast<ptrdiff_t>(write_index));
    return;
  }

  std::copy(source.begin() + static_cast<ptrdiff_t>(right_index), source.begin() + static_cast<ptrdiff_t>(right.end),
            destination.begin() + static_cast<ptrdiff_t>(write_index));
}

void MelnikIRadixSortIntTBB::MergeSortedRanges(std::vector<int> &data, std::vector<int> &buffer,
                                               const std::vector<Range> &ranges) {
  std::vector<int> *source = &data;
  std::vector<int> *destination = &buffer;
  std::vector<Range> current_ranges = ranges;

  while (current_ranges.size() > 1U) {
    const std::size_t merged_count = (current_ranges.size() + 1U) / 2U;
    std::vector<Range> next_ranges(merged_count);
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, merged_count),
                      [&](const tbb::blocked_range<std::size_t> &pair_block) {
      for (std::size_t pair_index = pair_block.begin(); pair_index < pair_block.end(); ++pair_index) {
        const std::size_t left_pos = pair_index * 2U;
        const Range left = current_ranges[left_pos];

        if (left_pos + 1U >= current_ranges.size()) {
          std::copy(source->begin() + static_cast<ptrdiff_t>(left.begin),
                    source->begin() + static_cast<ptrdiff_t>(left.end),
                    destination->begin() + static_cast<ptrdiff_t>(left.begin));
          next_ranges[pair_index] = left;
          continue;
        }

        const Range right = current_ranges[left_pos + 1U];
        MergeRanges(*source, *destination, left, right, left.begin);
        next_ranges[pair_index] = Range{.begin = left.begin, .end = right.end};
      }
    });

    current_ranges = std::move(next_ranges);
    std::swap(source, destination);
  }

  if (source != &data) {
    data.swap(*source);
  }
}

}  // namespace melnik_i_radix_sort_int
