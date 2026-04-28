#include "rozenberg_a_quicksort_simple_merge/tbb/include/ops_tbb.hpp"

#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/global_control.h>
#include <oneapi/tbb/parallel_for.h>
#include <tbb/tbb.h>

#include <stack>
#include <utility>
#include <vector>

#include "rozenberg_a_quicksort_simple_merge/common/include/common.hpp"

namespace rozenberg_a_quicksort_simple_merge {

RozenbergAQuicksortSimpleMergeTBB::RozenbergAQuicksortSimpleMergeTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());

  InType empty;
  GetInput().swap(empty);

  for (const auto &elem : in) {
    GetInput().push_back(elem);
  }

  GetOutput().clear();
}

bool RozenbergAQuicksortSimpleMergeTBB::ValidationImpl() {
  return (!(GetInput().empty())) && (GetOutput().empty());
}

bool RozenbergAQuicksortSimpleMergeTBB::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return GetOutput().size() == GetInput().size();
}

std::pair<int, int> RozenbergAQuicksortSimpleMergeTBB::Partition(InType &data, int left, int right) {
  const int pivot = data[left + ((right - left) / 2)];
  int i = left;
  int j = right;

  while (i <= j) {
    while (data[i] < pivot) {
      i++;
    }
    while (data[j] > pivot) {
      j--;
    }

    if (i <= j) {
      std::swap(data[i], data[j]);
      i++;
      j--;
    }
  }
  return {i, j};
}

void RozenbergAQuicksortSimpleMergeTBB::PushSubarrays(std::stack<std::pair<int, int>> &stack, int left, int right,
                                                      int i, int j) {
  if (j - left > right - i) {
    if (left < j) {
      stack.emplace(left, j);
    }
    if (i < right) {
      stack.emplace(i, right);
    }
  } else {
    if (i < right) {
      stack.emplace(i, right);
    }
    if (left < j) {
      stack.emplace(left, j);
    }
  }
}

void RozenbergAQuicksortSimpleMergeTBB::Quicksort(InType &data, int low, int high) {
  if (low >= high) {
    return;
  }

  std::stack<std::pair<int, int>> stack;

  stack.emplace(low, high);

  while (!stack.empty()) {
    const auto [left, right] = stack.top();
    stack.pop();

    if (left < right) {
      const auto [i, j] = Partition(data, left, right);
      PushSubarrays(stack, left, right, i, j);
    }
  }
}

void RozenbergAQuicksortSimpleMergeTBB::Merge(InType &data, int left, int mid, int right) {
  std::vector<int> temp(right - left + 1);
  int i = left;
  int j = mid + 1;
  int k = 0;

  while (i <= mid && j <= right) {
    if (data[i] <= data[j]) {
      temp[k++] = data[i++];
    } else {
      temp[k++] = data[j++];
    }
  }

  while (i <= mid) {
    temp[k++] = data[i++];
  }
  while (j <= right) {
    temp[k++] = data[j++];
  }

  for (int idx = 0; idx < k; ++idx) {
    data[left + idx] = temp[idx];
  }
}

bool RozenbergAQuicksortSimpleMergeTBB::RunImpl() {
  InType data = GetInput();
  int n = static_cast<int>(data.size());
  int num_threads = tbb::info::default_concurrency();

  tbb::global_control control(tbb::global_control::max_allowed_parallelism, num_threads);

  if (n < num_threads * 2) {
    Quicksort(data, 0, n - 1);
    GetOutput() = data;
    return true;
  }

  //  Create chunk borders container
  std::vector<int> borders(num_threads + 1);
  int chunk_size = n / num_threads;
  for (int i = 0; i < num_threads; i++) {
    borders[i] = i * chunk_size;
  }
  borders[num_threads] = n;

  //  Sort local chunks
  tbb::parallel_for(tbb::blocked_range<int>(0, num_threads), [&](const tbb::blocked_range<int> &range) {
    for (int i = range.begin(); i != range.end(); i++) {
      Quicksort(data, borders[i], borders[i + 1] - 1);
    }
  });

  //  Merge sorted chunks
  for (int i = 1; i < num_threads; i++) {
    Merge(data, 0, borders[i] - 1, borders[i + 1] - 1);
  }

  GetOutput() = data;
  return true;
}

bool RozenbergAQuicksortSimpleMergeTBB::PostProcessingImpl() {
  return true;
}

}  // namespace rozenberg_a_quicksort_simple_merge
