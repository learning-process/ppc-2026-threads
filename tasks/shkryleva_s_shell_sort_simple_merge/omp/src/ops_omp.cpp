#include "shkryleva_s_shell_sort_simple_merge/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <vector>

#include "shkryleva_s_shell_sort_simple_merge/common/include/common.hpp"

namespace shkryleva_s_shell_sort_simple_merge {

ShkrylevaSShellMergeOMP::ShkrylevaSShellMergeOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool ShkrylevaSShellMergeOMP::ValidationImpl() {
  return true;
}

bool ShkrylevaSShellMergeOMP::PreProcessingImpl() {
  input_data_ = GetInput();
  output_data_.clear();
  return true;
}

void ShkrylevaSShellMergeOMP::ShellSort(int left, int right, std::vector<int> &arr) {
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

void ShkrylevaSShellMergeOMP::Merge(int left, int mid, int right, std::vector<int> &arr, std::vector<int> &buffer) {
  int i = left;
  int j = mid + 1;
  int k = 0;

  int merge_size = right - left + 1;

  if (buffer.size() < static_cast<size_t>(merge_size)) {
    buffer.resize(static_cast<size_t>(merge_size));
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

bool ShkrylevaSShellMergeOMP::RunImpl() {
  if (input_data_.empty()) {
    output_data_.clear();
    return true;
  }

  std::vector<int> arr = input_data_;

  int array_size = static_cast<int>(arr.size());
  int num_threads = omp_get_max_threads();

  int sub_arr_size = (array_size + num_threads - 1) / num_threads;

#pragma omp parallel default(none) shared(arr, array_size, num_threads, sub_arr_size)
  {
    std::vector<int> buffer;

#pragma omp for schedule(dynamic)
    for (int i = 0; i < num_threads; ++i) {
      int left = i * sub_arr_size;
      int right = std::min(left + sub_arr_size - 1, array_size - 1);

      if (left < right) {
        ShellSort(left, right, arr);
      }
    }

    for (int size = sub_arr_size; size < array_size; size *= 2) {
#pragma omp for schedule(dynamic)
      for (int left = 0; left < array_size; left += 2 * size) {
        int mid = std::min(left + size - 1, array_size - 1);
        int right = std::min(left + (2 * size) - 1, array_size - 1);

        if (mid < right) {
          Merge(left, mid, right, arr, buffer);
        }
      }
    }
  }

  output_data_ = arr;

  return true;
}

bool ShkrylevaSShellMergeOMP::PostProcessingImpl() {
  GetOutput() = output_data_;
  return true;
}

}  // namespace shkryleva_s_shell_sort_simple_merge
