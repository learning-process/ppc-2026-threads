#include "shkryleva_s_shell_sort_simple_merge/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstddef>
#include <thread>
#include <vector>

#include "shkryleva_s_shell_sort_simple_merge/common/include/common.hpp"

namespace shkryleva_s_shell_sort_simple_merge {

ShkrylevaSShellMergeSTL::ShkrylevaSShellMergeSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool ShkrylevaSShellMergeSTL::ValidationImpl() {
  return true;
}

bool ShkrylevaSShellMergeSTL::PreProcessingImpl() {
  input_data_ = GetInput();
  output_data_.clear();
  return true;
}

void ShkrylevaSShellMergeSTL::ShellSort(int left, int right, std::vector<int> &arr) {
  int sub_array_size = right - left + 1;
  int gap = 1;
  while (gap <= sub_array_size / 3) {
    gap = (gap * 3) + 1;
  }
  for (; gap > 0; gap /= 3) {
    for (int i = left + gap; i <= right; ++i) {
      int temp = arr[i];
      int j = i;
      while (j >= left + gap && arr[j - gap] > temp) {
        arr[j] = arr[j - gap];
        j -= gap;
      }
      arr[j] = temp;
    }
  }
}

void ShkrylevaSShellMergeSTL::Merge(int left, int mid, int right, std::vector<int> &arr, std::vector<int> &buffer) {
  int i = left;
  int j = mid + 1;
  int k = 0;
  int merge_size = right - left + 1;
  if (static_cast<std::size_t>(merge_size) > buffer.size()) {
    buffer.resize(static_cast<std::size_t>(merge_size));
  }
  while (i <= mid || j <= right) {
    if (i > mid) {
      buffer[k++] = arr[j++];
    } else if (j > right) {
      buffer[k++] = arr[i++];
    } else {
      buffer[k++] = (arr[i] <= arr[j]) ? arr[i++] : arr[j++];
    }
  }
  for (int idx = 0; idx < k; ++idx) {
    arr[left + idx] = buffer[idx];
  }
}

void ShkrylevaSShellMergeSTL::SortSegments(std::vector<int> &arr, int num_threads, int sub_arr_size) {
  std::vector<std::thread> threads;
  for (int i = 0; i < num_threads; ++i) {
    int left = i * sub_arr_size;
    int right = std::min(left + sub_arr_size - 1, static_cast<int>(arr.size()) - 1);
    if (left < right) {
      threads.emplace_back([&arr, left, right] { ShellSort(left, right, arr); });
    }
  }
  for (auto &t : threads) {
    if (t.joinable()) {
      t.join();
    }
  }
}

void ShkrylevaSShellMergeSTL::HierarchicalMerge(std::vector<int> &arr, int num_threads, int sub_arr_size) {
  while (num_threads > 1) {
    int new_num_threads = (num_threads + 1) / 2;
    std::vector<std::thread> threads;
    for (int i = 0; i < new_num_threads; ++i) {
      int left = i * 2 * sub_arr_size;
      int mid = std::min(left + sub_arr_size - 1, static_cast<int>(arr.size()) - 1);
      int right = std::min(left + (2 * sub_arr_size) - 1, static_cast<int>(arr.size()) - 1);
      if (mid < right) {
        threads.emplace_back([&arr, left, mid, right] {
          std::vector<int> local_buffer;
          Merge(left, mid, right, arr, local_buffer);
        });
      }
    }
    for (auto &t : threads) {
      if (t.joinable()) {
        t.join();
      }
    }
    sub_arr_size *= 2;
    num_threads = new_num_threads;
  }
}

bool ShkrylevaSShellMergeSTL::RunImpl() {
  if (input_data_.empty()) {
    output_data_.clear();
    return true;
  }

  std::vector<int> arr = input_data_;
  const int array_size = static_cast<int>(arr.size());

  unsigned int hardware_threads = std::thread::hardware_concurrency();
  int num_threads = (hardware_threads > 0) ? static_cast<int>(hardware_threads) : 1;
  num_threads = std::min(num_threads, array_size);

  int sub_arr_size = (array_size + num_threads - 1) / num_threads;

  SortSegments(arr, num_threads, sub_arr_size);
  HierarchicalMerge(arr, num_threads, sub_arr_size);

  output_data_ = arr;
  return true;
}

bool ShkrylevaSShellMergeSTL::PostProcessingImpl() {
  GetOutput() = output_data_;
  return true;
}

}  // namespace shkryleva_s_shell_sort_simple_merge
