#include "rozenberg_a_quicksort_simple_merge/stl/include/ops_stl.hpp"

#include <stack>
#include <thread>
#include <utility>
#include <vector>

#include "rozenberg_a_quicksort_simple_merge/common/include/common.hpp"

namespace rozenberg_a_quicksort_simple_merge {

RozenbergAQuicksortSimpleMergeSTL::RozenbergAQuicksortSimpleMergeSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());

  InType empty;
  GetInput().swap(empty);

  for (const auto &elem : in) {
    GetInput().push_back(elem);
  }

  GetOutput().clear();
}

bool RozenbergAQuicksortSimpleMergeSTL::ValidationImpl() {
  return (!(GetInput().empty())) && (GetOutput().empty());
}

bool RozenbergAQuicksortSimpleMergeSTL::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return GetOutput().size() == GetInput().size();
}

std::pair<int, int> RozenbergAQuicksortSimpleMergeSTL::Partition(InType &data, int left, int right) {
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

void RozenbergAQuicksortSimpleMergeSTL::PushSubarrays(std::stack<std::pair<int, int>> &stack, int left, int right,
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

void RozenbergAQuicksortSimpleMergeSTL::Quicksort(InType &data, int low, int high) {
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

void RozenbergAQuicksortSimpleMergeSTL::Merge(InType &data, int left, int mid, int right) {
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

bool RozenbergAQuicksortSimpleMergeSTL::RunImpl() {
  InType data = GetInput();
  int n = static_cast<int>(data.size());

  int num_threads = static_cast<int>(std::thread::hardware_concurrency());
  if (num_threads == 0) {
    num_threads = 2;
  }

  if (n < num_threads * 2) {
    Quicksort(data, 0, n - 1);
    GetOutput() = data;
    return true;
  }

  //  Create chunk borders container
  std::vector<std::thread> threads;
  threads.reserve(num_threads);
  std::vector<int> borders(num_threads + 1);
  int chunk_size = n / num_threads;
  for (int i = 0; i < num_threads; i++) {
    borders[i] = i * chunk_size;
  }
  borders[num_threads] = n;

  //  Sort local chunks
  for (int i = 0; i < num_threads; i++) {
    threads.emplace_back([&data, &borders, i]() {
      rozenberg_a_quicksort_simple_merge::RozenbergAQuicksortSimpleMergeSTL::Quicksort(data, borders[i],
                                                                                       borders[i + 1] - 1);
    });
  }

  for (auto &t : threads) {
    if (t.joinable()) {
      t.join();
    }
  }

  //  Merge sorted chunks
  for (int i = 1; i < num_threads; i++) {
    Merge(data, 0, borders[i] - 1, borders[i + 1] - 1);
  }

  GetOutput() = data;
  return true;
}

bool RozenbergAQuicksortSimpleMergeSTL::PostProcessingImpl() {
  return true;
}

}  // namespace rozenberg_a_quicksort_simple_merge
