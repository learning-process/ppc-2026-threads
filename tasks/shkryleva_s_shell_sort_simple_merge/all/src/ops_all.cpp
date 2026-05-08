#include "shkryleva_s_shell_sort_simple_merge/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <cstddef>
#include <ranges>
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

// ========== MPI-вспомогательные функции ==========
std::vector<int> MergeTwoSorted(const std::vector<int> &left, const std::vector<int> &right) {
  std::vector<int> result;
  result.reserve(left.size() + right.size());
  size_t i = 0;
  size_t j = 0;
  while (i < left.size() && j < right.size()) {
    if (left[i] <= right[j]) {
      result.push_back(left[i++]);
    } else {
      result.push_back(right[j++]);
    }
  }
  while (i < left.size()) {
    result.push_back(left[i++]);
  }
  while (j < right.size()) {
    result.push_back(right[j++]);
  }
  return result;
}

void ExchangeAndMerge(int partner, std::vector<int> &merged_data) {
  size_t my_size = merged_data.size();
  size_t partner_size = 0;
  MPI_Sendrecv(&my_size, 1, MPI_UNSIGNED_LONG, partner, 0, &partner_size, 1, MPI_UNSIGNED_LONG, partner, 0,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  std::vector<int> partner_data(partner_size);
  MPI_Sendrecv(merged_data.data(), static_cast<int>(my_size), MPI_INT, partner, 1, partner_data.data(),
               static_cast<int>(partner_size), MPI_INT, partner, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  merged_data = MergeTwoSorted(merged_data, partner_data);
}

void ParallelHypercubeMerge(std::vector<int> &merged_data, int mpi_rank, int mpi_size) {
  int step = 1;
  while (step < mpi_size) {
    int partner = mpi_rank ^ step;
    if (partner < mpi_size) {
      ExchangeAndMerge(partner, merged_data);
    }
    step <<= 1;
  }
}

void BcastSortedVector(std::vector<int> &data, int mpi_rank) {
  size_t n = data.size();
  MPI_Bcast(&n, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
  if (mpi_rank != 0) {
    data.resize(n);
  }
  if (n > 0) {
    MPI_Bcast(data.data(), static_cast<int>(n), MPI_INT, 0, MPI_COMM_WORLD);
  }
}

void ComputeChunkParams(size_t total_size, int mpi_size, std::vector<size_t> &chunk_sizes,
                        std::vector<size_t> &offsets) {
  chunk_sizes.assign(mpi_size, 0);
  offsets.assign(mpi_size, 0);
  size_t base = total_size / static_cast<size_t>(mpi_size);
  size_t remainder = total_size % static_cast<size_t>(mpi_size);
  for (int i = 0; i < mpi_size; ++i) {
    chunk_sizes[static_cast<size_t>(i)] = base + (static_cast<size_t>(i) < static_cast<size_t>(remainder) ? 1U : 0U);
    offsets[static_cast<size_t>(i)] =
        (i == 0) ? 0 : offsets[static_cast<size_t>(i - 1)] + chunk_sizes[static_cast<size_t>(i - 1)];
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

  std::vector<size_t> chunk_sizes;
  std::vector<size_t> offsets;
  ComputeChunkParams(total_size, mpi_size, chunk_sizes, offsets);

  std::vector<int> local_data(chunk_sizes[static_cast<size_t>(mpi_rank)]);
  ScatterData(data, local_data, chunk_sizes, offsets);

  SortVectorParallel(local_data);

  if (mpi_size == 1) {
    data = std::move(local_data);
    return std::ranges::is_sorted(data);
  }

  std::vector<int> merged_data = std::move(local_data);
  ParallelHypercubeMerge(merged_data, mpi_rank, mpi_size);

  if (mpi_rank == 0) {
    data = std::move(merged_data);
  } else {
    data.clear();
  }

  BcastSortedVector(data, mpi_rank);

  MPI_Barrier(MPI_COMM_WORLD);
  return std::ranges::is_sorted(data);
}

bool ShkrylevaSShellMergeALL::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace shkryleva_s_shell_sort_simple_merge
