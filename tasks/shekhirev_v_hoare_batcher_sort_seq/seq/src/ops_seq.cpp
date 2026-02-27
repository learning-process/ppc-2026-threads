#include "shekhirev_v_hoare_batcher_sort_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <climits>
#include <cstddef>
#include <stack>
#include <tuple>
#include <utility>
#include <vector>

#include "shekhirev_v_hoare_batcher_sort_seq/common/include/common.hpp"

namespace shekhirev_v_hoare_batcher_sort_seq {

ShekhirevHoareBatcherSortSEQ::ShekhirevHoareBatcherSortSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool ShekhirevHoareBatcherSortSEQ::ValidationImpl() {
  return true;
}

bool ShekhirevHoareBatcherSortSEQ::PreProcessingImpl() {
  const auto &in = GetInput();
  size_t original_size = in.size();

  size_t p2 = 1;
  while (p2 < original_size) {
    p2 *= 2;
  }

  GetOutput() = in;
  GetOutput().resize(p2, INT_MAX);
  return true;
}

void ShekhirevHoareBatcherSortSEQ::HoareSort(std::vector<int> &arr, int left, int right) {
  if (left >= right) {
    return;
  }

  std::stack<std::pair<int, int>> tasks;
  tasks.emplace(left, right);

  while (!tasks.empty()) {
    auto [l, r] = tasks.top();
    tasks.pop();

    if (l >= r) {
      continue;
    }

    int pivot = arr[l + ((r - l) / 2)];
    int i = l;
    int j = r;

    while (i <= j) {
      while (arr[i] < pivot) {
        i++;
      }
      while (arr[j] > pivot) {
        j--;
      }
      if (i <= j) {
        std::swap(arr[i], arr[j]);
        i++;
        j--;
      }
    }

    if (i < r) {
      tasks.emplace(i, r);
    }
    if (l < j) {
      tasks.emplace(l, j);
    }
  }
}

void ShekhirevHoareBatcherSortSEQ::BatcherMerge(std::vector<int> &arr, int left, int right, int step) {
  struct MergeTask {
    int l;
    int r;
    int s;
    bool process;
  };

  std::stack<MergeTask> tasks;
  tasks.push({left, right, step, false});

  while (!tasks.empty()) {
    MergeTask task = tasks.top();
    tasks.pop();

    int n = task.r - task.l + 1;
    if (n <= task.s) {
      continue;
    }

    if (task.process) {
      for (int i = task.l + task.s; i + task.s <= task.r; i += task.s * 2) {
        if (arr[i] > arr[i + task.s]) {
          std::swap(arr[i], arr[i + task.s]);
        }
      }
    } else {
      tasks.push({task.l, task.r, task.s, true});

      tasks.push({task.l + task.s, task.r, task.s * 2, false});
      tasks.push({task.l, task.r, task.s * 2, false});
    }
  }
}

bool ShekhirevHoareBatcherSortSEQ::RunImpl() {
  auto &data = GetOutput();
  if (data.size() <= 1) {
    return true;
  }

  int n = static_cast<int>(data.size());
  int mid = n / 2;

  HoareSort(data, 0, mid - 1);
  HoareSort(data, mid, n - 1);

  BatcherMerge(data, 0, n - 1, 1);

  return true;
}

bool ShekhirevHoareBatcherSortSEQ::PostProcessingImpl() {
  size_t original_size = GetInput().size();
  GetOutput().resize(original_size);
  return true;
}

}  // namespace shekhirev_v_hoare_batcher_sort_seq
