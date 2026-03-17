#include "kondakov_v_shell_sort/omp/include/ops_omp.hpp"

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

#include "kondakov_v_shell_sort/common/include/common.hpp"
#include "util/include/util.hpp"

namespace kondakov_v_shell_sort {

namespace {

void ShellSortRange(std::vector<int> &data, size_t begin, size_t end) {
  const size_t n = (end > begin) ? (end - begin) : 0;
  for (size_t gap = n / 2; gap > 0; gap /= 2) {
    for (size_t i = begin + gap; i < end; ++i) {
      int value = data[i];
      size_t j = i;
      while (j >= begin + gap && data[j - gap] > value) {
        data[j] = data[j - gap];
        j -= gap;
      }
      data[j] = value;
    }
  }
}

void MergeRanges(const std::vector<int> &src, size_t left_begin, size_t left_end, size_t right_begin, size_t right_end,
                 std::vector<int> &dst, size_t dst_begin) {
  size_t i = left_begin;
  size_t j = right_begin;
  size_t k = dst_begin;

  while (i < left_end && j < right_end) {
    if (src[i] <= src[j]) {
      dst[k++] = src[i++];
    } else {
      dst[k++] = src[j++];
    }
  }

  while (i < left_end) {
    dst[k++] = src[i++];
  }
  while (j < right_end) {
    dst[k++] = src[j++];
  }
}

std::vector<std::pair<size_t, size_t>> MergeSegmentsOnce(const std::vector<int> &src, std::vector<int> &dst,
                                                         const std::vector<std::pair<size_t, size_t>> &segments,
                                                         int requested_threads) {
  const size_t pair_count = segments.size() / 2;
#pragma omp parallel for default(none) shared(segments, src, dst, pair_count) num_threads(requested_threads) \
    schedule(static)
  for (size_t pair = 0; pair < pair_count; ++pair) {
    const size_t base = pair * 2;
    const auto [left_begin, left_end] = segments[base];
    const auto [right_begin, right_end] = segments[base + 1];
    MergeRanges(src, left_begin, left_end, right_begin, right_end, dst, left_begin);
  }

  if ((segments.size() % 2) == 1) {
    const auto [tail_begin, tail_end] = segments.back();
    std::copy(src.begin() + static_cast<std::ptrdiff_t>(tail_begin),
              src.begin() + static_cast<std::ptrdiff_t>(tail_end),
              dst.begin() + static_cast<std::ptrdiff_t>(tail_begin));
  }

  std::vector<std::pair<size_t, size_t>> next_segments;
  next_segments.reserve(pair_count + (segments.size() % 2));
  for (size_t pair = 0; pair < pair_count; ++pair) {
    const size_t base = pair * 2;
    next_segments.emplace_back(segments[base].first, segments[base + 1].second);
  }
  if ((segments.size() % 2) == 1) {
    next_segments.push_back(segments.back());
  }
  return next_segments;
}

void MergeAllSegments(std::vector<int> &data, std::vector<int> &buffer, std::vector<std::pair<size_t, size_t>> segments,
                      int requested_threads) {
  bool sorted_in_data = true;
  while (segments.size() > 1) {
    const auto &src = sorted_in_data ? data : buffer;
    auto &dst = sorted_in_data ? buffer : data;
    segments = MergeSegmentsOnce(src, dst, segments, requested_threads);
    sorted_in_data = !sorted_in_data;
  }

  if (!sorted_in_data) {
    data.swap(buffer);
  }
}

}  // namespace

KondakovVShellSortOMP::KondakovVShellSortOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = in;
}

bool KondakovVShellSortOMP::ValidationImpl() {
  return !GetInput().empty();
}

bool KondakovVShellSortOMP::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool KondakovVShellSortOMP::RunImpl() {
  std::vector<int> &data = GetOutput();
  const size_t n = data.size();
  if (n <= 1) {
    return true;
  }

  const int requested_threads = std::max(1, ppc::util::GetNumThreads());
  const size_t parts = std::min(n, static_cast<size_t>(requested_threads));

  std::vector<size_t> offsets(parts + 1);
  for (size_t i = 0; i <= parts; ++i) {
    offsets[i] = (i * n) / parts;
  }

#pragma omp parallel for default(none) shared(data, offsets, parts) num_threads(requested_threads) schedule(static)
  for (size_t part = 0; part < parts; ++part) {
    const size_t begin = offsets[part];
    const size_t end = offsets[part + 1];
    ShellSortRange(data, begin, end);
  }

  std::vector<int> buffer(n);
  std::vector<std::pair<size_t, size_t>> segments;
  segments.reserve(parts);
  for (size_t part = 0; part < parts; ++part) {
    segments.emplace_back(offsets[part], offsets[part + 1]);
  }

  MergeAllSegments(data, buffer, std::move(segments), requested_threads);

  return std::ranges::is_sorted(data);
}

bool KondakovVShellSortOMP::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace kondakov_v_shell_sort
