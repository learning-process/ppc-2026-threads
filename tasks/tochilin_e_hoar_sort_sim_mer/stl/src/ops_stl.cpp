#include "tochilin_e_hoar_sort_sim_mer/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <thread>
#include <utility>
#include <vector>

#include "tochilin_e_hoar_sort_sim_mer/common/include/common.hpp"
#include "util/include/util.hpp"

namespace tochilin_e_hoar_sort_sim_mer {

namespace {

constexpr int kMinSequentialCutoff = 2048;

int ResolveConcurrency() {
  return std::max(1, ppc::util::GetNumThreads());
}

}  // namespace

TochilinEHoarSortSimMerSTL::TochilinEHoarSortSimMerSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool TochilinEHoarSortSimMerSTL::ValidationImpl() {
  return !GetInput().empty();
}

bool TochilinEHoarSortSimMerSTL::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

std::pair<int, int> TochilinEHoarSortSimMerSTL::Partition(std::vector<int> &arr, int l, int r) {
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

void TochilinEHoarSortSimMerSTL::QuickSortSequential(std::vector<int> &arr, int low, int high) {
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

bool TochilinEHoarSortSimMerSTL::ProcessRange(std::vector<int> &arr, std::pair<int, int> range, int serial_cutoff,
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

void TochilinEHoarSortSimMerSTL::QuickSortParallel(std::vector<int> &arr, int low, int high, int serial_cutoff,
                                                   int worker_count) {
  if (low >= high) {
    return;
  }

  std::vector<std::pair<int, int>> current_ranges;
  current_ranges.emplace_back(low, high);

  while (!current_ranges.empty()) {
    const int active_workers = std::max(1, std::min(worker_count, static_cast<int>(current_ranges.size())));
    std::vector<std::vector<std::pair<int, int>>> next_ranges_local(static_cast<std::size_t>(active_workers));
    std::vector<std::thread> workers;
    workers.reserve(static_cast<std::size_t>(active_workers));

    for (int worker_idx = 0; worker_idx < active_workers; ++worker_idx) {
      workers.emplace_back([&, worker_idx] {
        auto &local_next = next_ranges_local[static_cast<std::size_t>(worker_idx)];
        for (auto idx = static_cast<std::size_t>(worker_idx); idx < current_ranges.size();
             idx += static_cast<std::size_t>(active_workers)) {
          ProcessRange(arr, current_ranges[idx], serial_cutoff, local_next);
        }
      });
    }

    for (auto &worker : workers) {
      worker.join();
    }

    std::vector<std::pair<int, int>> next_ranges;
    for (auto &local_ranges : next_ranges_local) {
      next_ranges.insert(next_ranges.end(), local_ranges.begin(), local_ranges.end());
    }
    current_ranges = std::move(next_ranges);
  }
}

int TochilinEHoarSortSimMerSTL::ResolveSerialCutoff(std::size_t size) {
  const int concurrency = ResolveConcurrency();
  const std::size_t per_worker = size / static_cast<std::size_t>(concurrency * 4);
  return std::max(kMinSequentialCutoff, static_cast<int>(per_worker));
}

std::vector<int> TochilinEHoarSortSimMerSTL::MergeSortedVectors(const std::vector<int> &a, const std::vector<int> &b) {
  std::vector<int> result;
  result.reserve(a.size() + b.size());
  std::ranges::merge(a, b, std::back_inserter(result));
  return result;
}

bool TochilinEHoarSortSimMerSTL::RunImpl() {
  auto &data = GetOutput();
  if (data.empty()) {
    return false;
  }

  const auto mid = static_cast<std::vector<int>::difference_type>(data.size() / 2);
  const int concurrency = ResolveConcurrency();
  const int serial_cutoff = ResolveSerialCutoff(data.size());

  std::vector<int> left(data.begin(), data.begin() + mid);
  std::vector<int> right(data.begin() + mid, data.end());

  if (concurrency == 1) {
    QuickSortSequential(left, 0, static_cast<int>(left.size()) - 1);
    QuickSortSequential(right, 0, static_cast<int>(right.size()) - 1);
  } else {
    const int left_workers = std::max(1, concurrency / 2);
    const int right_workers = std::max(1, concurrency - left_workers);

    std::vector<std::thread> workers;
    workers.reserve(2);

    workers.emplace_back(
        [&]() { QuickSortParallel(left, 0, static_cast<int>(left.size()) - 1, serial_cutoff, left_workers); });
    workers.emplace_back(
        [&]() { QuickSortParallel(right, 0, static_cast<int>(right.size()) - 1, serial_cutoff, right_workers); });

    for (auto &worker : workers) {
      worker.join();
    }
  }

  data = MergeSortedVectors(left, right);

  return true;
}

bool TochilinEHoarSortSimMerSTL::PostProcessingImpl() {
  return std::ranges::is_sorted(GetOutput());
}

}  // namespace tochilin_e_hoar_sort_sim_mer
