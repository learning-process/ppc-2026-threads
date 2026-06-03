#include "karpich_i_bitwise_batcher/tbb/include/ops_tbb.hpp"

#include <algorithm>
#include <random>
#include <utility>
#include <vector>

#include "karpich_i_bitwise_batcher/common/include/common.hpp"
#include "oneapi/tbb/parallel_for.h"
#include "oneapi/tbb/parallel_invoke.h"

namespace karpich_i_bitwise_batcher {

namespace {

void RadixSortPositive(std::vector<int> &arr) {
  int n = static_cast<int>(arr.size());
  if (n <= 1) {
    return;
  }

  int max_val = *std::ranges::max_element(arr);
  if (max_val == 0) {
    return;
  }

  std::vector<int> buffer(n);

  for (int shift = 0; shift < 32 && (max_val >> shift) > 0; shift += 8) {
    std::vector<int> count(256, 0);
    for (int i = 0; i < n; i++) {
      count[(arr[i] >> shift) & 0xFF]++;
    }
    for (int i = 1; i < 256; i++) {
      count[i] += count[i - 1];
    }
    for (int i = n - 1; i >= 0; i--) {
      buffer[--count[(arr[i] >> shift) & 0xFF]] = arr[i];
    }
    arr = buffer;
  }
}

void RadixSort(std::vector<int> &arr) {
  int n = static_cast<int>(arr.size());
  if (n <= 1) {
    return;
  }

  std::vector<int> negative;
  std::vector<int> positive;
  for (int i = 0; i < n; i++) {
    if (arr[i] < 0) {
      negative.push_back(-arr[i]);
    } else {
      positive.push_back(arr[i]);
    }
  }

  RadixSortPositive(positive);
  RadixSortPositive(negative);

  int idx = 0;
  for (int i = static_cast<int>(negative.size()) - 1; i >= 0; i--) {
    arr[idx++] = -negative[i];
  }
  for (int x : positive) {
    arr[idx++] = x;
  }
}

struct MergeTask {
  int lo;
  int hi;
  int r;
};

std::vector<std::vector<std::pair<int, int>>> BuildMergeNetwork(int lo, int hi) {
  std::vector<std::vector<std::pair<int, int>>> levels;
  std::vector<MergeTask> current = {{lo, hi, 1}};

  while (!current.empty()) {
    std::vector<MergeTask> next;
    std::vector<std::pair<int, int>> comps;

    for (const auto &[tlo, thi, tr] : current) {
      int step = tr * 2;
      if (step < thi - tlo) {
        next.push_back({tlo, thi, step});
        next.push_back({tlo + tr, thi, step});
        for (int i = tlo + tr; i + tr <= thi; i += step) {
          comps.emplace_back(i, i + tr);
        }
      } else if (tlo + tr <= thi) {
        comps.emplace_back(tlo, tlo + tr);
      }
    }

    levels.push_back(std::move(comps));
    current = std::move(next);
  }

  return levels;
}

void BatcherMergeParallel(std::vector<int> &arr, int lo, int hi) {
  auto levels = BuildMergeNetwork(lo, hi);
  for (int lvl = static_cast<int>(levels.size()) - 1; lvl >= 0; lvl--) {
    const auto &level = levels[lvl];
    int level_size = static_cast<int>(level.size());
    tbb::parallel_for(0, level_size, [&arr, &level](int idx) {
      auto [a, b] = level[idx];
      if (arr[a] > arr[b]) {
        std::swap(arr[a], arr[b]);
      }
    });
  }
}

}  // namespace

KarpichIBitwiseBatcherTBB::KarpichIBitwiseBatcherTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool KarpichIBitwiseBatcherTBB::ValidationImpl() {
  return GetInput() > 0;
}

bool KarpichIBitwiseBatcherTBB::PreProcessingImpl() {
  int n = GetInput();
  data_.resize(n);
  std::mt19937 gen(static_cast<unsigned int>(n));
  std::uniform_int_distribution<int> dist(-1000, 1000);
  for (int i = 0; i < n; i++) {
    data_[i] = dist(gen);
  }
  return true;
}

bool KarpichIBitwiseBatcherTBB::RunImpl() {
  int n = static_cast<int>(data_.size());
  if (n <= 1) {
    return true;
  }

  int padded = 1;
  while (padded < n) {
    padded *= 2;
  }

  int max_elem = *std::ranges::max_element(data_);
  data_.resize(padded, max_elem);

  int half = padded / 2;
  std::vector<int> left(data_.begin(), data_.begin() + half);
  std::vector<int> right(data_.begin() + half, data_.end());

  tbb::parallel_invoke([&left]() { RadixSort(left); }, [&right]() { RadixSort(right); });

  std::ranges::copy(left, data_.begin());
  std::ranges::copy(right, data_.begin() + half);

  BatcherMergeParallel(data_, 0, padded - 1);

  data_.resize(n);
  return true;
}

bool KarpichIBitwiseBatcherTBB::PostProcessingImpl() {
  for (int i = 1; std::cmp_less(i, data_.size()); i++) {
    if (data_[i] < data_[i - 1]) {
      return false;
    }
  }
  GetOutput() = GetInput();
  return true;
}

}  // namespace karpich_i_bitwise_batcher
