#include "../include/radix_sort_omp.hpp"

#include <omp.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <utility>
#include <vector>

#include "util/include/util.hpp"

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
}  // namespace

void RadixSortOMP::SortRange(std::vector<double> &arr, std::size_t left, std::size_t right) {
  if (right <= left) {
    return;
  }

  std::size_t n = right - left;
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
      auto b = static_cast<std::size_t>((data[i] >> (byte * 8)) & kByteMask);
      ++count.at(b);
    }

    std::size_t sum = 0;
    for (std::size_t i = 0; i < kRadix; ++i) {
      std::size_t tmp = count.at(i);
      count.at(i) = sum;
      sum += tmp;
    }

    for (std::size_t i = 0; i < n; ++i) {
      auto b = static_cast<std::size_t>((data[i] >> (byte * 8)) & kByteMask);
      buffer[count.at(b)++] = data[i];
    }

    data.swap(buffer);
  }

  for (std::size_t i = 0; i < n; ++i) {
    std::uint64_t bits = FromSortable(data[i]);
    std::memcpy(&arr[left + i], &bits, sizeof(double));
  }
}

std::vector<double> RadixSortOMP::Merge(const std::vector<double> &a, const std::vector<double> &b) {
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

void RadixSortOMP::Sort(std::vector<double> &arr) {
  if (arr.empty()) {
    return;
  }

  int num_threads = ppc::util::GetNumThreads();
  if (num_threads <= 0) {
    num_threads = 1;
  }

  if (static_cast<std::size_t>(num_threads) > arr.size()) {
    num_threads = static_cast<int>(arr.size());
  }

  std::vector<std::pair<std::size_t, std::size_t>> ranges(num_threads);

  std::size_t base = arr.size() / static_cast<std::size_t>(num_threads);
  std::size_t rem = arr.size() % static_cast<std::size_t>(num_threads);

  std::size_t begin = 0;
  for (std::size_t i = 0; i < ranges.size(); ++i) {
    std::size_t block_size = base + (i < rem ? 1 : 0);
    ranges[i] = {begin, begin + block_size};
    begin += block_size;
  }

#pragma omp parallel default(none) shared(arr, ranges, num_threads) num_threads(ppc::util::GetNumThreads())
  {
#pragma omp for schedule(static)
    for (int i = 0; i < num_threads; ++i) {
      SortRange(arr, ranges[i].first, ranges[i].second);
    }
  }

  std::vector<std::vector<double>> parts(num_threads);
  for (int i = 0; i < num_threads; ++i) {
    parts[i] = std::vector<double>(arr.begin() + static_cast<std::ptrdiff_t>(ranges[i].first),
                                   arr.begin() + static_cast<std::ptrdiff_t>(ranges[i].second));
  }

  while (parts.size() > 1) {
    std::size_t pair_count = parts.size() / 2;
    const int pair_count_int = static_cast<int>(pair_count);
    std::vector<std::vector<double>> next((parts.size() + 1) / 2);

#pragma omp parallel default(none) shared(parts, next, pair_count_int) num_threads(ppc::util::GetNumThreads())
    {
#pragma omp for schedule(static)
      for (int i = 0; i < pair_count_int; ++i) {
        const auto idx = static_cast<std::size_t>(i);
        next[idx] = Merge(parts[2 * idx], parts[(2 * idx) + 1]);
      }
    }

    if (parts.size() % 2 != 0) {
      next.back() = std::move(parts.back());
    }

    parts = std::move(next);
  }

  arr = std::move(parts[0]);
}
