#include "tochilin_e_hoar_sort_sim_mer/tbb/include/ops_tbb.hpp"

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <utility>
#include <vector>

#include "oneapi/tbb/enumerable_thread_specific.h"
#include "oneapi/tbb/global_control.h"
#include "oneapi/tbb/parallel_for.h"
#include "oneapi/tbb/parallel_invoke.h"
#include "tochilin_e_hoar_sort_sim_mer/common/include/common.hpp"
#include "util/include/util.hpp"

namespace tochilin_e_hoar_sort_sim_mer {

namespace {

constexpr int kMinSequentialCutoff = 2048;

int ResolveConcurrency() {
  return std::max(1, ppc::util::GetNumThreads());
}

}  // namespace

TochilinEHoarSortSimMerTBB::TochilinEHoarSortSimMerTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool TochilinEHoarSortSimMerTBB::ValidationImpl() {
  return !GetInput().empty();
}

bool TochilinEHoarSortSimMerTBB::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

std::pair<int, int> TochilinEHoarSortSimMerTBB::Partition(std::vector<int> &arr, int l, int r) {
  int i = l;
  int j = r;
  const int pivot = arr[(l + r) / 2];

  while (i <= j) {
    while (arr[i] < pivot) {
      ++i;
    }
    while (arr[j] > pivot) {
      --j;
    }
    if (i <= j) {
      std::swap(arr[i], arr[j]);
      ++i;
      --j;
    }
  }

  return {i, j};
}

void TochilinEHoarSortSimMerTBB::QuickSortSequential(std::vector<int> &arr, int low, int high) {
  if (low >= high) {
    return;
  }

  std::vector<std::pair<int, int>> stack;
  stack.emplace_back(low, high);

  while (!stack.empty()) {
    const auto [l, r] = stack.back();
    stack.pop_back();

    if (l >= r) {
      continue;
    }

    const auto [i, j] = Partition(arr, l, r);

    if (l < j) {
      stack.emplace_back(l, j);
    }
    if (i < r) {
      stack.emplace_back(i, r);
    }
  }
}

bool TochilinEHoarSortSimMerTBB::ProcessRange(std::vector<int> &arr, std::pair<int, int> range, int serial_cutoff,
                                              std::vector<std::pair<int, int>> &next_ranges) {
  const auto [left, right] = range;
  if (left >= right) {
    return false;
  }

  const int range_size = right - left + 1;
  if (range_size <= serial_cutoff) {
    QuickSortSequential(arr, left, right);
    return false;
  }

  const auto [i, j] = Partition(arr, left, right);
  if (left < j) {
    next_ranges.emplace_back(left, j);
  }
  if (i < right) {
    next_ranges.emplace_back(i, right);
  }

  return true;
}

void TochilinEHoarSortSimMerTBB::QuickSortParallel(std::vector<int> &arr, int low, int high, int serial_cutoff) {
  if (low >= high) {
    return;
  }

  std::vector<std::pair<int, int>> current_ranges;
  current_ranges.emplace_back(low, high);

  while (!current_ranges.empty()) {
    tbb::enumerable_thread_specific<std::vector<std::pair<int, int>>> next_ranges_local;

    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, current_ranges.size()),
                      [&](const tbb::blocked_range<std::size_t> &range) {
      auto &local_next = next_ranges_local.local();
      for (std::size_t idx = range.begin(); idx != range.end(); ++idx) {
        ProcessRange(arr, current_ranges[idx], serial_cutoff, local_next);
      }
    });

    std::vector<std::pair<int, int>> next_ranges;
    for (auto &local_ranges : next_ranges_local) {
      next_ranges.insert(next_ranges.end(), local_ranges.begin(), local_ranges.end());
    }
    current_ranges = std::move(next_ranges);
  }
}

int TochilinEHoarSortSimMerTBB::ResolveSerialCutoff(std::size_t size) {
  const int concurrency = ResolveConcurrency();
  const std::size_t per_worker = size / static_cast<std::size_t>(concurrency * 4);
  return std::max(kMinSequentialCutoff, static_cast<int>(per_worker));
}

std::vector<int> TochilinEHoarSortSimMerTBB::MergeSortedVectors(const std::vector<int> &a, const std::vector<int> &b) {
  std::vector<int> result;
  result.reserve(a.size() + b.size());
  std::ranges::merge(a, b, std::back_inserter(result));
  return result;
}

bool TochilinEHoarSortSimMerTBB::RunImpl() {
  auto &data = GetOutput();
  if (data.empty()) {
    return false;
  }

  const tbb::global_control control(tbb::global_control::max_allowed_parallelism,
                                    static_cast<std::size_t>(ResolveConcurrency()));
  const auto mid = static_cast<std::vector<int>::difference_type>(data.size() / 2);
  const int serial_cutoff = ResolveSerialCutoff(data.size());

  std::vector<int> left(data.begin(), data.begin() + mid);
  std::vector<int> right(data.begin() + mid, data.end());

  tbb::parallel_invoke([&]() { QuickSortParallel(left, 0, static_cast<int>(left.size()) - 1, serial_cutoff); },
                       [&]() { QuickSortParallel(right, 0, static_cast<int>(right.size()) - 1, serial_cutoff); });

  data = MergeSortedVectors(left, right);

  return true;
}

bool TochilinEHoarSortSimMerTBB::PostProcessingImpl() {
  return std::ranges::is_sorted(GetOutput());
}

}  // namespace tochilin_e_hoar_sort_sim_mer
