#include "shkryleva_s_shell_sort_simple_merge/tbb/include/ops_tbb.hpp"

#include <oneapi/tbb/info.h>
#include <oneapi/tbb/task_arena.h>
#include <oneapi/tbb/task_group.h>

#include <algorithm>
#include <cstddef>
#include <thread>
#include <utility>
#include <vector>

#include "shkryleva_s_shell_sort_simple_merge/common/include/common.hpp"
namespace shkryleva_s_shell_sort_simple_merge {

ShkrylevaSShellMergeTBB::ShkrylevaSShellMergeTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool ShkrylevaSShellMergeTBB::ValidationImpl() {
  return true;
}

bool ShkrylevaSShellMergeTBB::PreProcessingImpl() {
  input_data_ = GetInput();
  output_data_.clear();
  return true;
}

void ShkrylevaSShellMergeTBB::ShellSort(int left, int right, std::vector<int> &arr) {
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

void ShkrylevaSShellMergeTBB::Merge(int left, int mid, int right, std::vector<int> &arr, std::vector<int> &buffer) {
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

bool ShkrylevaSShellMergeTBB::RunImpl() {
  if (input_data_.empty()) {
    output_data_.clear();
    return true;
  }

  std::vector<int> arr = input_data_;
  const int array_size = static_cast<int>(arr.size());

  // Определение максимального количества потоков
  // Используем TBB, если доступно, иначе std::thread
  int max_threads = tbb::info::default_concurrency();
  if (max_threads <= 0) {
    max_threads = static_cast<int>(std::thread::hardware_concurrency());
    if (max_threads <= 0) {
      max_threads = 4;
    }
  }

  int threads = std::min(max_threads, array_size);
  const int sub_arr_size = (array_size + threads - 1) / threads;

  // Разбиение на сегменты
  std::vector<std::pair<int, int>> segments;
  segments.reserve(threads);
  for (int idx = 0; idx < threads; ++idx) {
    int left = idx * sub_arr_size;
    int right = std::min(left + sub_arr_size - 1, array_size - 1);
    segments.emplace_back(left, right);
  }

  // Параллельная сортировка каждого сегмента
  tbb::task_arena arena(threads);
  arena.execute([&] {
    tbb::task_group tg;
    for (const auto &seg : segments) {
      int l = seg.first;
      int r = seg.second;
      tg.run([&arr, l, r] { ShellSort(l, r, arr); });
    }
    tg.wait();
  });

  // Последовательное слияние отсортированных сегментов
  std::vector<int> buffer;
  int current_end = segments.front().second;
  for (std::size_t i = 1; i < segments.size(); ++i) {
    int next_end = segments[i].second;
    Merge(0, current_end, next_end, arr, buffer);
    current_end = next_end;
  }

  output_data_ = arr;
  return true;
}

bool ShkrylevaSShellMergeTBB::PostProcessingImpl() {
  GetOutput() = output_data_;
  return true;
}

}  // namespace shkryleva_s_shell_sort_simple_merge
