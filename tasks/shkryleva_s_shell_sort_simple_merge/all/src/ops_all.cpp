#include "shkryleva_s_shell_sort_simple_merge/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <cstddef>
#include <thread>
#include <utility>
#include <vector>

#include "shkryleva_s_shell_sort_simple_merge/common/include/common.hpp"

namespace shkryleva_s_shell_sort_simple_merge {

namespace {

// ========== Локальная параллельная сортировка (STL-версия) ==========
void ShellSort(int left, int right, std::vector<int> &arr) {
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

void Merge(int left, int mid, int right, std::vector<int> &arr, std::vector<int> &buffer) {
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

void SortSegments(std::vector<int> &arr, int num_threads, int sub_arr_size) {
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

void HierarchicalMerge(std::vector<int> &arr, int num_threads, int sub_arr_size) {
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

void SortVectorParallel(std::vector<int> &arr) {
  if (arr.size() < 2) {
    return;
  }

  const int array_size = static_cast<int>(arr.size());
  unsigned int hardware_threads = std::thread::hardware_concurrency();
  int num_threads = (hardware_threads > 0) ? static_cast<int>(hardware_threads) : 1;
  num_threads = std::min(num_threads, array_size);

  int sub_arr_size = (array_size + num_threads - 1) / num_threads;

  SortSegments(arr, num_threads, sub_arr_size);
  HierarchicalMerge(arr, num_threads, sub_arr_size);
}

void ComputeChunkParams(size_t total_size, int mpi_size, std::vector<size_t> &chunk_sizes,
                        std::vector<size_t> &offsets) {
  const auto mpi_size_usz = static_cast<size_t>(mpi_size);
  chunk_sizes.assign(mpi_size_usz, 0);
  offsets.assign(mpi_size_usz, 0);
  const size_t base = total_size / mpi_size_usz;
  const size_t remainder = total_size % mpi_size_usz;
  for (size_t i = 0; i < mpi_size_usz; ++i) {
    chunk_sizes[i] = base + (i < remainder ? 1U : 0U);
    offsets[i] = (i == 0) ? 0 : offsets[i - 1] + chunk_sizes[i - 1];
  }
}

void ScatterData(const std::vector<int> &global_data, std::vector<int> &local_data,
                 const std::vector<size_t> &chunk_sizes, const std::vector<size_t> &offsets) {
  int mpi_size = static_cast<int>(chunk_sizes.size());
  std::vector<int> send_counts(mpi_size);
  std::vector<int> send_displs(mpi_size);
  for (int i = 0; i < mpi_size; ++i) {
    send_counts[i] = static_cast<int>(chunk_sizes[static_cast<size_t>(i)]);
    send_displs[i] = static_cast<int>(offsets[static_cast<size_t>(i)]);
  }
  MPI_Scatterv(global_data.data(), send_counts.data(), send_displs.data(), MPI_INT, local_data.data(),
               static_cast<int>(local_data.size()), MPI_INT, 0, MPI_COMM_WORLD);
}

}  // namespace

// ========== Реализация класса ShkrylevaSShellMergeALL ==========
ShkrylevaSShellMergeALL::ShkrylevaSShellMergeALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = in;
}

bool ShkrylevaSShellMergeALL::ValidationImpl() {
  return !GetInput().empty();
}

bool ShkrylevaSShellMergeALL::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool ShkrylevaSShellMergeALL::RunImpl() {
  int mpi_rank = 0;
  int mpi_size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  std::vector<int> &data = GetOutput();
  const size_t total_size = data.size();
  if (total_size <= 1) {
    return true;
  }

  // Разбиение на чанки
  std::vector<size_t> chunk_sizes;
  std::vector<size_t> offsets;
  ComputeChunkParams(total_size, mpi_size, chunk_sizes, offsets);

  std::vector<int> local_data(chunk_sizes[mpi_rank]);
  ScatterData(data, local_data, chunk_sizes, offsets);

  // Локальная параллельная сортировка (STL-потоки)
  SortVectorParallel(local_data);

  // Сбор всех отсортированных кусков на процессе 0
  std::vector<int> all_data;
  if (mpi_rank == 0) {
    all_data.resize(total_size);
    std::vector<int> recv_counts(mpi_size);
    for (int i = 0; i < mpi_size; ++i) {
      recv_counts[i] = static_cast<int>(chunk_sizes[i]);
    }
    std::vector<int> displs(mpi_size, 0);
    for (int i = 1; i < mpi_size; ++i) {
      displs[i] = displs[i - 1] + recv_counts[i - 1];
    }
    MPI_Gatherv(local_data.data(), static_cast<int>(local_data.size()), MPI_INT, all_data.data(), recv_counts.data(),
                displs.data(), MPI_INT, 0, MPI_COMM_WORLD);
  } else {
    MPI_Gatherv(local_data.data(), static_cast<int>(local_data.size()), MPI_INT, nullptr, nullptr, nullptr, MPI_INT, 0,
                MPI_COMM_WORLD);
  }

  // На процессе 0 – финальная сортировка (можно std::sort)
  if (mpi_rank == 0) {
    std::ranges::sort(all_data);
    data = std::move(all_data);
  }

  // Рассылка результата всем процессам
  if (mpi_rank == 0) {
    size_t sz = data.size();
    MPI_Bcast(&sz, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(data.data(), static_cast<int>(sz), MPI_INT, 0, MPI_COMM_WORLD);
  } else {
    size_t sz = 0;
    MPI_Bcast(&sz, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    data.resize(sz);
    MPI_Bcast(data.data(), static_cast<int>(sz), MPI_INT, 0, MPI_COMM_WORLD);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  return std::ranges::is_sorted(data);
}

bool ShkrylevaSShellMergeALL::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace shkryleva_s_shell_sort_simple_merge
