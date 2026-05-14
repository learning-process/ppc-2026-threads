#include "shkryleva_s_shell_sort_simple_merge/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <vector>

#include "shkryleva_s_shell_sort_simple_merge/common/include/common.hpp"

namespace shkryleva_s_shell_sort_simple_merge {

// ==================== Реализация вспомогательных функций ====================

void ShkrylevaSShellMergeALL::ShellSort(std::vector<int> &arr, int left, int right) {
  int sub_size = right - left + 1;
  int gap = 1;
  while (gap <= sub_size / 3) {
    gap = gap * 3 + 1;
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
  if (num_threads > n) {
    num_threads = n;
  }

  // Разбиение на блоки
  int block_size = (n + num_threads - 1) / num_threads;
  std::vector<std::vector<int>> local_buffers(num_threads);  // буферы для слияния

  // Параллельная сортировка блоков
#pragma omp parallel for schedule(static) num_threads(num_threads)
  for (int t = 0; t < num_threads; ++t) {
    int left = t * block_size;
    int right = std::min(left + block_size - 1, n - 1);
    if (left < right) {
      ShellSort(arr, left, right);
    }
  }

  // Иерархическое слияние блоков
  int step = block_size;
  int active = num_threads;
  while (active > 1) {
    int new_active = (active + 1) / 2;
#pragma omp parallel for schedule(static) num_threads(new_active)
    for (int t = 0; t < new_active; ++t) {
      int left = t * 2 * step;
      int mid = std::min(left + step - 1, n - 1);
      int right = std::min(left + 2 * step - 1, n - 1);
      if (mid < right) {
        Merge(arr, left, mid, right, local_buffers[t]);
      }
    }
    step *= 2;
    active = new_active;
  }
}

std::vector<int> ShkrylevaSShellMergeALL::SimpleMerge(const std::vector<int> &a, const std::vector<int> &b) {
  std::vector<int> res;
  res.reserve(a.size() + b.size());
  size_t i = 0, j = 0;
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

// ==================== Реализация методов класса ====================

ShkrylevaSShellMergeALL::ShkrylevaSShellMergeALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<int>();  // как в успешной задаче
}

bool ShkrylevaSShellMergeALL::ValidationImpl() {
  return !GetInput().empty();
}

bool ShkrylevaSShellMergeALL::PreProcessingImpl() {
  return true;  // ничего не делаем
}

bool ShkrylevaSShellMergeALL::RunImpl() {
  int rank = 0, size_comm = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size_comm);

  int total_size = static_cast<int>(GetInput().size());
  MPI_Bcast(&total_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (total_size == 0) {
    return true;
  }

  // Подготовка параметров для Scatterv
  std::vector<int> send_counts(size_comm), offsets(size_comm);
  int chunk = total_size / size_comm;
  int rem = total_size % size_comm;
  for (int i = 0; i < size_comm; ++i) {
    send_counts[i] = chunk + (i < rem ? 1 : 0);
    offsets[i] = (i == 0) ? 0 : offsets[i - 1] + send_counts[i - 1];
  }

  std::vector<int> local_data(send_counts[rank]);
  const int *in_ptr = (rank == 0) ? GetInput().data() : nullptr;
  MPI_Scatterv(const_cast<int *>(in_ptr), send_counts.data(), offsets.data(), MPI_INT, local_data.data(),
               send_counts[rank], MPI_INT, 0, MPI_COMM_WORLD);

  // Локальная сортировка с использованием OpenMP
  ParallelShellSort(local_data);

  // Сбор и слияние на процессе 0
  if (rank == 0) {
    std::vector<int> final = local_data;
    for (int i = 1; i < size_comm; ++i) {
      if (send_counts[i] == 0) {
        continue;
      }
      std::vector<int> recv_buf(send_counts[i]);
      MPI_Recv(recv_buf.data(), send_counts[i], MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      final = SimpleMerge(final, recv_buf);
    }
    GetOutput() = std::move(final);
  } else {
    if (!local_data.empty()) {
      MPI_Send(local_data.data(), send_counts[rank], MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
  }

  // Рассылка отсортированного массива всем процессам
  GetOutput().resize(static_cast<std::size_t>(total_size));
  MPI_Bcast(GetOutput().data(), total_size, MPI_INT, 0, MPI_COMM_WORLD);

  return true;
}

bool ShkrylevaSShellMergeALL::PostProcessingImpl() {
  return true;
}

}  // namespace shkryleva_s_shell_sort_simple_merge
