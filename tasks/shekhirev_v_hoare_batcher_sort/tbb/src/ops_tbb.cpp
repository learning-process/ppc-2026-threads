#include "shekhirev_v_hoare_batcher_sort/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <algorithm>
#include <climits>
#include <utility>
#include <vector>

#include "shekhirev_v_hoare_batcher_sort/common/include/common.hpp"
#include "util/include/util.hpp"

namespace shekhirev_v_hoare_batcher_sort {

namespace {

void HoareSort(std::vector<int> &arr, int low, int high) {
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
      while (arr[++i] < pivot) {
      }
      while (arr[--j] > pivot) {
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

void ProcessBatcherStep(std::vector<int> &a, int step_p, int step_k, int start, int total_pairs) {
  tbb::parallel_for(tbb::blocked_range<int>(0, total_pairs), [&](const tbb::blocked_range<int> &r) {
    for (int step = r.begin(); step != r.end(); ++step) {
      int i = step % step_k;
      int b = step / step_k;
      int j = start + (b * (step_k * 2));

      int idx1 = i + j;
      int idx2 = i + j + step_k;

      if ((idx1 / (step_p * 2)) == (idx2 / (step_p * 2))) {
        if (a[idx1] > a[idx2]) {
          std::swap(a[idx1], a[idx2]);
        }
      }
    }
  });
}

void BatcherMerge(std::vector<int> &a, int n_pow2, int chunk_size) {
  for (int step_p = chunk_size; step_p < n_pow2; step_p *= 2) {
    for (int step_k = step_p; step_k >= 1; step_k /= 2) {
      int start = step_k;
      int num_blocks = (n_pow2 / (2 * step_k)) - 1;

      if (step_k == step_p) {
        start = 0;
        num_blocks = n_pow2 / (2 * step_k);
      }

      int total_pairs = num_blocks * step_k;
      ProcessBatcherStep(a, step_p, step_k, start, total_pairs);
    }
  }
}

}  // namespace

ShekhirevHoareBatcherSortTBB::ShekhirevHoareBatcherSortTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool ShekhirevHoareBatcherSortTBB::ValidationImpl() {
  return true;
}

bool ShekhirevHoareBatcherSortTBB::PreProcessingImpl() {
  input_ = GetInput();
  return true;
}

bool ShekhirevHoareBatcherSortTBB::RunImpl() {
  int n = static_cast<int>(input_.size());
  if (n <= 1) {
    res_.assign(input_.begin(), input_.end());
    return true;
  }

  int n_pow2 = 1;
  while (n_pow2 < n) {
    n_pow2 *= 2;
  }

  std::vector<int> a(n_pow2, INT_MAX);
  for (int i = 0; i < n; i++) {
    a[i] = input_[i];
  }

  int num_threads = ppc::util::GetNumThreads();
  tbb::global_control global_limit(tbb::global_control::max_allowed_parallelism, num_threads);

  int p_threads = 1;
  while ((p_threads * 2) <= num_threads && (p_threads * 2) <= n_pow2) {
    p_threads *= 2;
  }

  int chunk_size = n_pow2 / p_threads;

  tbb::parallel_for(tbb::blocked_range<int>(0, p_threads), [&](const tbb::blocked_range<int> &r) {
    for (int i = r.begin(); i != r.end(); ++i) {
      HoareSort(a, i * chunk_size, ((i + 1) * chunk_size) - 1);
    }
  });

  BatcherMerge(a, n_pow2, chunk_size);

  res_.assign(a.begin(), a.begin() + n);
  return true;
}

bool ShekhirevHoareBatcherSortTBB::PostProcessingImpl() {
  GetOutput() = res_;
  return true;
}

}  // namespace shekhirev_v_hoare_batcher_sort
