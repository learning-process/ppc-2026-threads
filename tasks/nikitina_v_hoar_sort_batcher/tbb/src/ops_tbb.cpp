#include "nikitina_v_hoar_sort_batcher/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <algorithm>
#include <limits>
#include <utility>
#include <vector>

#include "nikitina_v_hoar_sort_batcher/common/include/common.hpp"

namespace nikitina_v_hoar_sort_batcher {

namespace {

void QuickSortHoare(std::vector<int> &arr, int low, int high) {
  if (low >= high) {
    return;
  }
  std::vector<std::pair<int, int>> stack;
  stack.emplace_back(low, high);

  while (!stack.empty()) {
    auto [l, h] = stack.back();
    stack.pop_back();

    if (l >= h) {
      continue;
    }

    int pivot = arr[l + ((h - l) / 2)];
    int i = l - 1;
    int j = h + 1;

    while (true) {
      i++;
      while (arr[i] < pivot) {
        i++;
      }

      j--;
      while (arr[j] > pivot) {
        j--;
      }

      if (i >= j) {
        break;
      }
      std::swap(arr[i], arr[j]);
    }

    stack.emplace_back(l, j);
    stack.emplace_back(j + 1, h);
  }
}

void CompareSplit(std::vector<int> &arr, int start1, int len1, int start2, int len2) {
  if (len1 == 0 || len2 == 0) {
    return;
  }

  std::vector<int> left_block(arr.begin() + start1, arr.begin() + start1 + len1);
  std::vector<int> right_block(arr.begin() + start2, arr.begin() + start2 + len2);

  int p1 = 0;
  int p2 = 0;
  int write1 = start1;
  int write2 = start2;

  for (int i = 0; i < len1 + len2; ++i) {
    int val = 0;
    if (p1 < len1 && (p2 == len2 || left_block[p1] <= right_block[p2])) {
      val = left_block[p1++];
    } else {
      val = right_block[p2++];
    }

    if (i < len1) {
      arr[write1++] = val;
    } else {
      arr[write2++] = val;
    }
  }
}

void BuildPairs(std::vector<std::pair<int, int>> &pairs, int num_threads, int step_p, int step_k) {
  for (int idx_j = step_k % step_p; idx_j + step_k < num_threads; idx_j += (step_k * 2)) {
    for (int idx_i = 0; idx_i < std::min(step_k, num_threads - idx_j - step_k); idx_i++) {
      if ((idx_j + idx_i) / (step_p * 2) == (idx_j + idx_i + step_k) / (step_p * 2)) {
        pairs.emplace_back(idx_j + idx_i, idx_j + idx_i + step_k);
      }
    }
  }
}

void BatcherMergePhase(std::vector<int> &output, const std::vector<int> &offsets, int num_threads) {
  for (int step_p = 1; step_p < num_threads; step_p *= 2) {
    for (int step_k = step_p; step_k > 0; step_k /= 2) {
      std::vector<std::pair<int, int>> pairs;
      BuildPairs(pairs, num_threads, step_p, step_k);

      int num_pairs = static_cast<int>(pairs.size());

      // TBB аналог: #pragma omp parallel for
      tbb::parallel_for(tbb::blocked_range<int>(0, num_pairs), [&](const tbb::blocked_range<int> &r) {
        for (int idx = r.begin(); idx != r.end(); ++idx) {
          int block_a = pairs[idx].first;
          int block_b = pairs[idx].second;
          CompareSplit(output, offsets[block_a], offsets[block_a + 1] - offsets[block_a], offsets[block_b],
                       offsets[block_b + 1] - offsets[block_b]);
        }
      });
    }
  }
}

}  // namespace

HoareSortBatcherTBB::HoareSortBatcherTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool HoareSortBatcherTBB::ValidationImpl() {
  return true;
}

bool HoareSortBatcherTBB::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool HoareSortBatcherTBB::RunImpl() {
  auto &output = GetOutput();

  int orig_n = static_cast<int>(output.size());
  if (orig_n <= 1) {
    return true;
  }

  int max_threads = tbb::this_task_arena::max_concurrency();
  int t = 1;
  while (t * 2 <= max_threads && t * 2 <= orig_n) {
    t *= 2;
  }

  if (t == 1) {
    QuickSortHoare(output, 0, orig_n - 1);
    return true;
  }

  int pad = (t - (orig_n % t)) % t;
  for (int i = 0; i < pad; ++i) {
    output.push_back(std::numeric_limits<int>::max());
  }

  int n = orig_n + pad;
  std::vector<int> offsets(t + 1, 0);
  int chunk = n / t;
  for (int i = 0; i <= t; ++i) {
    offsets[i] = i * chunk;
  }

  tbb::parallel_for(tbb::blocked_range<int>(0, t), [&](const tbb::blocked_range<int> &r) {
    for (int i = r.begin(); i != r.end(); ++i) {
      QuickSortHoare(output, offsets[i], offsets[i + 1] - 1);
    }
  });

  BatcherMergePhase(output, offsets, t);

  output.resize(orig_n);

  return true;
}

bool HoareSortBatcherTBB::PostProcessingImpl() {
  return true;
}

}  // namespace nikitina_v_hoar_sort_batcher
