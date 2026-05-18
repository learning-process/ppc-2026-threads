#include "shkryleva_s_shell_sort_simple_merge/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

#include "shkryleva_s_shell_sort_simple_merge/common/include/common.hpp"

namespace shkryleva_s_shell_sort_simple_merge {

// ========== вспомогательные функции для распределения и сборки ==========
namespace {

void ComputeDistribution(int total_size, int size_comm, std::vector<int> &send_counts, std::vector<int> &offsets) {
  int chunk = total_size / size_comm;
  int rem = total_size % size_comm;
  for (int i = 0; i < size_comm; ++i) {
    send_counts[i] = chunk + (i < rem ? 1 : 0);
    offsets[i] = (i == 0) ? 0 : offsets[i - 1] + send_counts[i - 1];
  }
}

void GatherAndMerge(int rank, int size_comm, const std::vector<int> &local_data, const std::vector<int> &send_counts,
                    std::vector<int> &final_result) {
  if (rank == 0) {
    final_result = local_data;
    for (int i = 1; i < size_comm; ++i) {
      if (send_counts[i] == 0) {
        continue;
      }
      std::vector<int> recv_buf(send_counts[i]);
      MPI_Recv(recv_buf.data(), send_counts[i], MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      final_result = ShkrylevaSShellMergeALL::SimpleMerge(final_result, recv_buf);
    }
  } else {
    if (!local_data.empty()) {
      MPI_Send(local_data.data(), send_counts[rank], MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
  }
}

}  // namespace

// ========== реализация методов сортировки ==========

void ShkrylevaSShellMergeALL::ShellSort(std::vector<int> &arr, int left, int right) {
  int sub_size = right - left + 1;
  int gap = 1;
  while (gap <= sub_size / 3) {
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

void ShkrylevaSShellMergeALL::Merge(std::vector<int> &arr, int left, int mid, int right, std::vector<int> &buffer) {
  int i = left;
  int j = mid + 1;
  int k = 0;
  int size = right - left + 1;
  if (static_cast<std::size_t>(size) > buffer.size()) {
    buffer.resize(static_cast<std::size_t>(size));
  }
  while (i <= mid && j <= right) {
    buffer[k++] = (arr[i] <= arr[j]) ? arr[i++] : arr[j++];
  }
  while (i <= mid) {
    buffer[k++] = arr[i++];
  }
  while (j <= right) {
    buffer[k++] = arr[j++];
  }
  for (int idx = 0; idx < k; ++idx) {
    arr[left + idx] = buffer[idx];
  }
}

void ShkrylevaSShellMergeALL::ParallelShellSort(std::vector<int> &arr) {
  int n = static_cast<int>(arr.size());
  if (n < 2) {
    return;
  }

  int num_threads = omp_get_max_threads();
  num_threads = std::min(num_threads, n);
  int block_size = (n + num_threads - 1) / num_threads;

// Параллельная сортировка блоков
#pragma omp parallel for schedule(static) num_threads(num_threads) default(none) shared(arr, n, block_size)
  for (int tid = 0; tid < num_threads; ++tid) {
    int left = tid * block_size;
    int right = std::min(left + block_size - 1, n - 1);
    if (left < right) {
      ShellSort(arr, left, right);
    }
  }

  // Последовательное слияние блоков (без OpenMP)
  std::vector<int> buffer;
  int step = block_size;
  while (step < n) {
    for (int left = 0; left < n; left += 2 * step) {
      int mid = std::min(left + step - 1, n - 1);
      int right = std::min(left + 2 * step - 1, n - 1);
      if (mid < right) {
        Merge(arr, left, mid, right, buffer);
      }
    }
    step *= 2;
  }
}

std::vector<int> ShkrylevaSShellMergeALL::SimpleMerge(const std::vector<int> &a, const std::vector<int> &b) {
  std::vector<int> res;
  res.reserve(a.size() + b.size());
  size_t i = 0;
  size_t j = 0;
  while (i < a.size() && j < b.size()) {
    if (a[i] <= b[j]) {
      res.push_back(a[i++]);
    } else {
      res.push_back(b[j++]);
    }
  }
  while (i < a.size()) {
    res.push_back(a[i++]);
  }
  while (j < b.size()) {
    res.push_back(b[j++]);
  }
  return res;
}

// ========== реализация методов класса ==========

ShkrylevaSShellMergeALL::ShkrylevaSShellMergeALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<int>();
}

bool ShkrylevaSShellMergeALL::ValidationImpl() {
  return !GetInput().empty();
}

bool ShkrylevaSShellMergeALL::PreProcessingImpl() {
  return true;
}

bool ShkrylevaSShellMergeALL::RunImpl() {
  int rank = 0;
  int size_comm = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size_comm);

  int total_size = static_cast<int>(GetInput().size());
  MPI_Bcast(&total_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (total_size == 0) {
    return true;
  }

  std::vector<int> send_counts(size_comm);
  std::vector<int> offsets(size_comm);
  ComputeDistribution(total_size, size_comm, send_counts, offsets);

  std::vector<int> local_data(send_counts[rank]);
  int *in_ptr = (rank == 0) ? GetInput().data() : nullptr;
  MPI_Scatterv(in_ptr, send_counts.data(), offsets.data(), MPI_INT, local_data.data(), send_counts[rank], MPI_INT, 0,
               MPI_COMM_WORLD);

  ParallelShellSort(local_data);

  std::vector<int> final_result;
  GatherAndMerge(rank, size_comm, local_data, send_counts, final_result);

  if (rank == 0) {
    GetOutput() = std::move(final_result);
  }

  GetOutput().resize(static_cast<std::size_t>(total_size));
  MPI_Bcast(GetOutput().data(), total_size, MPI_INT, 0, MPI_COMM_WORLD);

  return true;
}

bool ShkrylevaSShellMergeALL::PostProcessingImpl() {
  return true;
}

}  // namespace shkryleva_s_shell_sort_simple_merge
