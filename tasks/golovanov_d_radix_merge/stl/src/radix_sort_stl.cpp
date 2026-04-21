#include "../include/radix_sort_stl.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <thread>
#include <utility>
#include <vector>

#include "util/include/util.hpp"

namespace golovanov_d_radix_merge {

namespace {
constexpr int kBytes = 8;
constexpr std::size_t kRadix = 256;
constexpr std::uint64_t kSignMask = 1ULL << 63;
constexpr std::uint64_t kByteMask = 0xFFULL;

std::uint64_t ToSortable(std::uint64_t bits) {
  return ((bits & kSignMask) != 0U) ? ~bits : (bits ^ kSignMask);
}

std::uint64_t FromSortable(std::uint64_t bits) {
  return ((bits & kSignMask) != 0U) ? (bits ^ kSignMask) : ~bits;
}

int GetThreadCount(std::size_t size) {
  if (size == 0) {
    return 1;
  }

  int num_threads = ppc::util::GetNumThreads();
  if (num_threads <= 0) {
    num_threads = 1;
  }

  const auto max_threads = static_cast<int>(size);
  return std::min(num_threads, max_threads);
}

std::vector<std::pair<std::size_t, std::size_t>> BuildRanges(std::size_t size, int num_threads) {
  if (num_threads <= 0) {
    num_threads = 1;
  }

  const auto threads_count = static_cast<std::size_t>(num_threads);
  std::vector<std::pair<std::size_t, std::size_t>> ranges(threads_count);

  const std::size_t base = size / threads_count;
  const std::size_t rem = size % threads_count;

  std::size_t begin = 0;
  for (std::size_t i = 0; i < ranges.size(); ++i) {
    const std::size_t block_size = base + ((i < rem) ? 1U : 0U);
    ranges[i] = {begin, begin + block_size};
    begin += block_size;
  }

  return ranges;
}

void SortParts(std::vector<double> &arr, const std::vector<std::pair<std::size_t, std::size_t>> &ranges) {
  std::vector<std::thread> workers;
  workers.reserve(ranges.size());

  for (const auto &range : ranges) {
    workers.emplace_back([&arr, range]() { RadixSortSTL::SortRange(arr, range.first, range.second); });
  }

  for (auto &worker : workers) {
    worker.join();
  }
}

std::vector<std::vector<double>> CopyParts(const std::vector<double> &arr,
                                           const std::vector<std::pair<std::size_t, std::size_t>> &ranges) {
  std::vector<std::vector<double>> parts(ranges.size());

  for (std::size_t i = 0; i < ranges.size(); ++i) {
    parts[i] = std::vector<double>(arr.begin() + static_cast<std::ptrdiff_t>(ranges[i].first),
                                   arr.begin() + static_cast<std::ptrdiff_t>(ranges[i].second));
  }

  return parts;
}

std::vector<std::vector<double>> MergeStep(const std::vector<std::vector<double>> &parts) {
  const std::size_t pair_count = parts.size() / 2;
  std::vector<std::vector<double>> next((parts.size() + 1) / 2);
  std::vector<std::thread> workers;
  workers.reserve(pair_count);

  for (std::size_t i = 0; i < pair_count; ++i) {
    workers.emplace_back([&parts, &next, i]() { next[i] = RadixSortSTL::Merge(parts[2 * i], parts[(2 * i) + 1]); });
  }

  for (auto &worker : workers) {
    worker.join();
  }

  if ((parts.size() % 2U) != 0U) {
    next.back() = parts.back();
  }

  return next;
}

}  // namespace

void RadixSortSTL::SortRange(std::vector<double> &arr, std::size_t left, std::size_t right) {
  if (right <= left) {
    return;
  }

  const std::size_t n = right - left;
  std::vector<std::uint64_t> data(n);

  for (std::size_t i = 0; i < n; ++i) {
    std::uint64_t bits = 0;
    std::memcpy(&bits, &arr[left + i], sizeof(double));
    data[i] = ToSortable(bits);
  }

  std::vector<std::uint64_t> buffer(n);

  for (int byte = 0; byte < kBytes; ++byte) {
    std::array<std::size_t, kRadix> count{};

    for (std::size_t i = 0; i < n; ++i) {
      const auto bucket = static_cast<std::size_t>((data[i] >> (byte * 8)) & kByteMask);
      ++count.at(bucket);
    }

    std::size_t sum = 0;
    for (std::size_t i = 0; i < kRadix; ++i) {
      const std::size_t tmp = count.at(i);
      count.at(i) = sum;
      sum += tmp;
    }

    for (std::size_t i = 0; i < n; ++i) {
      const auto bucket = static_cast<std::size_t>((data[i] >> (byte * 8)) & kByteMask);
      const std::size_t pos = count.at(bucket);
      buffer[pos] = data[i];
      ++count.at(bucket);
    }

    data.swap(buffer);
  }

  for (std::size_t i = 0; i < n; ++i) {
    const std::uint64_t bits = FromSortable(data[i]);
    std::memcpy(&arr[left + i], &bits, sizeof(double));
  }
}

std::vector<double> RadixSortSTL::Merge(const std::vector<double> &a, const std::vector<double> &b) {
  std::vector<double> result;
  result.reserve(a.size() + b.size());

  std::size_t i = 0;
  std::size_t j = 0;

  while (i < a.size() && j < b.size()) {
    if (a[i] <= b[j]) {
      result.push_back(a[i]);
      ++i;
    } else {
      result.push_back(b[j]);
      ++j;
    }
  }

  while (i < a.size()) {
    result.push_back(a[i]);
    ++i;
  }

  while (j < b.size()) {
    result.push_back(b[j]);
    ++j;
  }

  return result;
}

void RadixSortSTL::Sort(std::vector<double> &arr) {
  if (arr.empty()) {
    return;
  }

  const int num_threads = GetThreadCount(arr.size());
  const auto ranges = BuildRanges(arr.size(), num_threads);

  SortParts(arr, ranges);

  std::vector<std::vector<double>> parts = CopyParts(arr, ranges);
  while (parts.size() > 1) {
    parts = MergeStep(parts);
  }

  arr = std::move(parts.front());
}

}  // namespace golovanov_d_radix_merge
