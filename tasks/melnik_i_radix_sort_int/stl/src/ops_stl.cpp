#include "melnik_i_radix_sort_int/stl/include/ops_stl.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <thread>
#include <utility>
#include <vector>

#include "melnik_i_radix_sort_int/common/include/common.hpp"
#include "util/include/util.hpp"

namespace melnik_i_radix_sort_int {

namespace {

std::vector<MelnikIRadixSortIntSTL::Range> BuildInitialRanges(std::size_t data_size, int num_ranges,
                                                              std::vector<int> *initial_buffer) {
  std::vector<MelnikIRadixSortIntSTL::Range> ranges;
  ranges.reserve(static_cast<std::size_t>(num_ranges));
  const std::size_t chunk_size =
      (data_size + static_cast<std::size_t>(num_ranges) - 1U) / static_cast<std::size_t>(num_ranges);

  for (int range_index = 0; range_index < num_ranges; ++range_index) {
    const std::size_t begin = static_cast<std::size_t>(range_index) * chunk_size;
    if (begin >= data_size) {
      break;
    }
    const std::size_t end = std::min(begin + chunk_size, data_size);
    ranges.push_back(MelnikIRadixSortIntSTL::Range{.begin = begin, .end = end, .buffer = initial_buffer});
  }

  return ranges;
}

}  // namespace

MelnikIRadixSortIntSTL::MelnikIRadixSortIntSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool MelnikIRadixSortIntSTL::ValidationImpl() {
  return !GetInput().empty();
}

bool MelnikIRadixSortIntSTL::PreProcessingImpl() {
  GetOutput() = GetInput();
  return !GetOutput().empty();
}

bool MelnikIRadixSortIntSTL::RunImpl() {
  if (GetOutput().empty()) {
    return false;
  }

  const std::size_t data_size = GetOutput().size();
  const int requested_threads = std::max(1, ppc::util::GetNumThreads());
  const int num_threads = std::min<int>(requested_threads, static_cast<int>(data_size));

  std::vector<int> buffer(data_size);
  auto &output = GetOutput();

  if (num_threads <= 1) {
    auto *final_buffer = RadixSortRange(output, buffer, 0, data_size);
    if (final_buffer != &output) {
      output.swap(*final_buffer);
    }
    return !output.empty();
  }

  std::vector<Range> ranges = BuildInitialRanges(data_size, num_threads, &output);
  const int active_ranges = static_cast<int>(ranges.size());

  std::vector<std::thread> workers;
  workers.reserve(static_cast<std::size_t>(active_ranges));

  for (int i = 0; i < active_ranges; ++i) {
    workers.emplace_back([&, i]() {
      Range &range = ranges[static_cast<std::size_t>(i)];
      if (range.begin < range.end) {
        range.buffer = RadixSortRange(output, buffer, range.begin, range.end);
      }
    });
  }

  for (auto &worker : workers) {
    worker.join();
  }

  MergeSortedRanges(output, buffer, ranges);
  return !GetOutput().empty();
}

bool MelnikIRadixSortIntSTL::PostProcessingImpl() {
  return std::ranges::is_sorted(GetOutput());
}

void MelnikIRadixSortIntSTL::CountingSortByByte(const std::vector<int> &source, std::vector<int> &destination,
                                                std::size_t begin, std::size_t end, std::int64_t exp,
                                                std::int64_t offset, const std::array<std::size_t, kBuckets> &count) {
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

std::vector<int> *MelnikIRadixSortIntSTL::RadixSortRange(std::vector<int> &data, std::vector<int> &buffer,
                                                         std::size_t begin, std::size_t end) {
  if (end - begin <= 1) {
    return &data;
  }

  const auto range_begin = data.begin() + static_cast<ptrdiff_t>(begin);
  const auto range_end = data.begin() + static_cast<ptrdiff_t>(end);
  const auto [min_it, max_it] = std::ranges::minmax_element(range_begin, range_end);
  const auto min_value = static_cast<std::int64_t>(*min_it);
  const auto max_value = static_cast<std::int64_t>(*max_it);
  const std::int64_t offset = (min_value < 0) ? -min_value : 0;
  const std::int64_t max_shifted_value = max_value + offset;

  std::vector<std::array<std::size_t, kBuckets>> all_counts;
  for (std::int64_t exp = 1; max_shifted_value / exp > 0; exp <<= kBitsPerPass) {
    all_counts.emplace_back();
    all_counts.back().fill(0);
  }

  for (std::size_t index = begin; index < end; ++index) {
    std::int64_t val = static_cast<std::int64_t>(data[index]) + offset;
    for (auto &count_table : all_counts) {
      count_table.at(static_cast<std::size_t>(val % static_cast<std::int64_t>(kBuckets)))++;
      val /= static_cast<std::int64_t>(kBuckets);
    }
  }

  std::vector<int> *source = &data;
  std::vector<int> *destination = &buffer;
  std::int64_t exp = 1;

  for (const auto &count_table : all_counts) {
    CountingSortByByte(*source, *destination, begin, end, exp, offset, count_table);
    std::swap(source, destination);
    exp <<= kBitsPerPass;
  }

  return source;
}

void MelnikIRadixSortIntSTL::MergeRanges(const std::vector<int> &left_source, const std::vector<int> &right_source,
                                         std::vector<int> &destination, Range left, Range right,
                                         std::size_t write_begin) {
  std::size_t left_index = left.begin;
  std::size_t right_index = right.begin;
  std::size_t write_index = write_begin;

  while (left_index < left.end && right_index < right.end) {
    if (left_source[left_index] <= right_source[right_index]) {
      destination[write_index] = left_source[left_index];
      ++left_index;
    } else {
      destination[write_index] = right_source[right_index];
      ++right_index;
    }
    ++write_index;
  }

  if (left_index < left.end) {
    std::copy(left_source.begin() + static_cast<ptrdiff_t>(left_index),
              left_source.begin() + static_cast<ptrdiff_t>(left.end),
              destination.begin() + static_cast<ptrdiff_t>(write_index));
  } else if (right_index < right.end) {
    std::copy(right_source.begin() + static_cast<ptrdiff_t>(right_index),
              right_source.begin() + static_cast<ptrdiff_t>(right.end),
              destination.begin() + static_cast<ptrdiff_t>(write_index));
  }
}

void MelnikIRadixSortIntSTL::EnsureRangeInBuffer(Range &range, std::vector<int> *target_buffer) {
  if (range.buffer == target_buffer) {
    return;
  }
  std::copy(range.buffer->begin() + static_cast<ptrdiff_t>(range.begin),
            range.buffer->begin() + static_cast<ptrdiff_t>(range.end),
            target_buffer->begin() + static_cast<ptrdiff_t>(range.begin));
  range.buffer = target_buffer;
}

MelnikIRadixSortIntSTL::Range MelnikIRadixSortIntSTL::ProcessMergePair(std::size_t pair_index,
                                                                       const std::vector<Range> &current_ranges,
                                                                       std::vector<int> *majority_source,
                                                                       std::vector<int> *level_dest) {
  const std::size_t left_pos = pair_index * 2U;
  Range left = current_ranges[left_pos];

  if (left_pos + 1U >= current_ranges.size()) {
    EnsureRangeInBuffer(left, majority_source);
    return left;
  }

  Range right = current_ranges[left_pos + 1U];
  EnsureRangeInBuffer(left, majority_source);
  EnsureRangeInBuffer(right, majority_source);
  MergeRanges(*majority_source, *majority_source, *level_dest, left, right, left.begin);
  return Range{.begin = left.begin, .end = right.end, .buffer = level_dest};
}

void MelnikIRadixSortIntSTL::MergeSortedRanges(std::vector<int> &data, std::vector<int> &buffer,
                                               std::vector<Range> &ranges) {
  if (ranges.empty()) {
    return;
  }

  std::vector<Range> current_ranges = ranges;

  while (current_ranges.size() > 1U) {
    const std::size_t merged_count = (current_ranges.size() + 1U) / 2U;
    std::vector<Range> next_ranges(merged_count);
    std::vector<int> *majority_source = current_ranges[0].buffer;
    std::vector<int> *level_dest = (majority_source == &data) ? &buffer : &data;

    std::vector<std::thread> workers;
    workers.reserve(merged_count);

    for (std::size_t pair_index = 0; pair_index < merged_count; ++pair_index) {
      workers.emplace_back([&, pair_index, majority_source, level_dest]() {
        next_ranges[pair_index] = ProcessMergePair(pair_index, current_ranges, majority_source, level_dest);
      });
    }

    for (auto &worker : workers) {
      worker.join();
    }

    current_ranges = std::move(next_ranges);
  }

  if (current_ranges[0].buffer != &data) {
    data.swap(*current_ranges[0].buffer);
  }
}

}  // namespace melnik_i_radix_sort_int
