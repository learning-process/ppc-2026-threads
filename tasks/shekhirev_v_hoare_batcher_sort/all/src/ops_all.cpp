#include "shekhirev_v_hoare_batcher_sort/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

#include "shekhirev_v_hoare_batcher_sort/common/include/common.hpp"

namespace shekhirev_v_hoare_batcher_sort {

namespace {

void SplitPartition(std::vector<int> &arr, int &l, int &r, int &i, int &j) {
  int pivot = arr[l + ((r - l) / 2)];
  i = l;
  j = r;
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
}

void ProcessPartition(std::vector<int> &arr, int &l, int &r, std::vector<std::pair<int, int>> &stack) {
  int i = 0;
  int j = 0;
  SplitPartition(arr, l, r, i, j);

  if (j - l < r - i) {
    if (i < r) {
      stack.emplace_back(i, r);
    }
    r = j;
  } else {
    if (l < j) {
      stack.emplace_back(l, j);
    }
    l = i;
  }
}

void OptimizedHoareSort(std::vector<int> &arr, int left, int right) {
  if (left >= right) {
    return;
  }
  std::vector<std::pair<int, int>> stack;
  stack.reserve(64);
  stack.emplace_back(left, right);

  while (!stack.empty()) {
    auto [l, r] = stack.back();
    stack.pop_back();
    while (l < r) {
      ProcessPartition(arr, l, r, stack);
    }
  }
}

void OMPLocalSort(std::vector<int> &arr) {
  int size = static_cast<int>(arr.size());
  if (size <= 1) {
    return;
  }

  int max_threads = omp_get_max_threads();
  if (max_threads <= 0) {
    max_threads = 1;
  }

  int num_threads = 1;
  while (num_threads * 2 <= max_threads && num_threads * 2 <= size) {
    num_threads *= 2;
  }

  if (num_threads <= 1) {
    OptimizedHoareSort(arr, 0, size - 1);
    return;
  }

#pragma omp parallel for default(none) shared(arr, num_threads, size)
  for (int i = 0; i < num_threads; ++i) {
    int chunk_size = size / num_threads;
    int start = i * chunk_size;
    int end = (i == num_threads - 1) ? size - 1 : (start + chunk_size - 1);
    OptimizedHoareSort(arr, start, end);
  }

  OptimizedHoareSort(arr, 0, size - 1);
}

void MpiCompareAndSwap(std::vector<int> &local_arr, int neighbor, bool keep_low_half) {
  int size = static_cast<int>(local_arr.size());
  std::vector<int> neighbor_arr(static_cast<std::size_t>(size));

  MPI_Sendrecv(local_arr.data(), size, MPI_INT, neighbor, 0, neighbor_arr.data(), size, MPI_INT, neighbor, 0,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  std::vector<int> merged_arr(static_cast<std::size_t>(size) * 2);
  std::ranges::merge(local_arr, neighbor_arr, merged_arr.begin());

  if (keep_low_half) {
    std::copy(merged_arr.begin(), merged_arr.begin() + size, local_arr.begin());
  } else {
    std::copy(merged_arr.begin() + size, merged_arr.end(), local_arr.begin());
  }
}

void BatcherExchange(std::vector<int> &local_data, int rank, int world_size, int p_step, int k_step) {
  for (int j_idx = k_step % p_step; j_idx + k_step < world_size; j_idx += (k_step * 2)) {
    int upper_bound = std::min(k_step, world_size - j_idx - k_step);
    for (int i_idx = 0; i_idx < upper_bound; ++i_idx) {
      int r1 = j_idx + i_idx;
      int r2 = j_idx + i_idx + k_step;
      if ((r1 / (p_step * 2)) == (r2 / (p_step * 2))) {
        if (rank == r1) {
          MpiCompareAndSwap(local_data, r2, true);
        } else if (rank == r2) {
          MpiCompareAndSwap(local_data, r1, false);
        }
      }
    }
  }
}

void BatcherMergeNetwork(std::vector<int> &local_data, int rank, int world_size) {
  for (int p_step = 1; p_step < world_size; p_step *= 2) {
    for (int k_step = p_step; k_step > 0; k_step /= 2) {
      BatcherExchange(local_data, rank, world_size, p_step, k_step);
    }
  }
}

}  // namespace

ShekhirevHoareBatcherSortALL::ShekhirevHoareBatcherSortALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool ShekhirevHoareBatcherSortALL::ValidationImpl() {
  return true;
}

bool ShekhirevHoareBatcherSortALL::PreProcessingImpl() {
  input_ = GetInput();
  return true;
}

bool ShekhirevHoareBatcherSortALL::RunImpl() {
  int rank = 0;
  int world_size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  int total_n = (rank == 0) ? static_cast<int>(input_.size()) : 0;
  MPI_Bcast(&total_n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (total_n <= 1) {
    if (rank == 0) {
      output_ = input_;
    }
    return true;
  }

  int chunk_size = (total_n + world_size - 1) / world_size;
  int padded_total_size = chunk_size * world_size;

  std::vector<int> local_data(static_cast<std::size_t>(chunk_size));
  std::vector<int> send_buffer;

  if (rank == 0) {
    send_buffer = input_;
    send_buffer.resize(static_cast<std::size_t>(padded_total_size), std::numeric_limits<int>::max());
  }

  MPI_Scatter(send_buffer.data(), chunk_size, MPI_INT, local_data.data(), chunk_size, MPI_INT, 0, MPI_COMM_WORLD);

  OMPLocalSort(local_data);

  BatcherMergeNetwork(local_data, rank, world_size);

  std::vector<int> gather_buffer;
  if (rank == 0) {
    gather_buffer.resize(static_cast<std::size_t>(padded_total_size));
  }
  MPI_Gather(local_data.data(), chunk_size, MPI_INT, gather_buffer.data(), chunk_size, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    gather_buffer.resize(total_n);
    output_ = std::move(gather_buffer);
  }

  return true;
}

bool ShekhirevHoareBatcherSortALL::PostProcessingImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    GetOutput() = output_;
  }
  return true;
}

}  // namespace shekhirev_v_hoare_batcher_sort
