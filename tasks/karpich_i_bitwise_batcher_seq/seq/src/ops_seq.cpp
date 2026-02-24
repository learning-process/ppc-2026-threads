#include "karpich_i_bitwise_batcher_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <random>
#include <utility>
#include <vector>

#include "karpich_i_bitwise_batcher_seq/common/include/common.hpp"

namespace karpich_i_bitwise_batcher_seq {

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

void ApplyComparatorNetwork(std::vector<int> &arr, const std::vector<std::vector<std::pair<int, int>>> &levels) {
  for (int lvl = static_cast<int>(levels.size()) - 1; lvl >= 0; lvl--) {
    for (const auto &[a, b] : levels[lvl]) {
      if (arr[a] > arr[b]) {
        std::swap(arr[a], arr[b]);
      }
    }
  }
}

void BatcherMerge(std::vector<int> &arr, int lo, int hi) {
  auto levels = BuildMergeNetwork(lo, hi);
  ApplyComparatorNetwork(arr, levels);
}

}  // namespace

KarpichIBitwiseBatcherSEQ::KarpichIBitwiseBatcherSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool KarpichIBitwiseBatcherSEQ::ValidationImpl() {
  return GetInput() > 0;
}

bool KarpichIBitwiseBatcherSEQ::PreProcessingImpl() {
  int n = GetInput();
  data_.resize(n);
  std::mt19937 gen(static_cast<unsigned int>(n));
  std::uniform_int_distribution<int> dist(-1000, 1000);
  for (int i = 0; i < n; i++) {
    data_[i] = dist(gen);
  }
  return true;
}

bool KarpichIBitwiseBatcherSEQ::RunImpl() {
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

  RadixSort(left);
  RadixSort(right);

  std::ranges::copy(left, data_.begin());
  std::ranges::copy(right, data_.begin() + half);

  BatcherMerge(data_, 0, padded - 1);

  data_.resize(n);
  return true;
}

bool KarpichIBitwiseBatcherSEQ::PostProcessingImpl() {
  for (int i = 1; std::cmp_less(i, data_.size()); i++) {
    if (data_[i] < data_[i - 1]) {
      return false;
    }
  }
  GetOutput() = GetInput();
  return true;
}

}  // namespace karpich_i_bitwise_batcher_seq
