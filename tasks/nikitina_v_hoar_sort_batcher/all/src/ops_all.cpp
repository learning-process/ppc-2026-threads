#include "nikitina_v_hoar_sort_batcher/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <cstddef>
#include <limits>
#include <thread>
#include <utility>
#include <vector>

#include "nikitina_v_hoar_sort_batcher/common/include/common.hpp"

namespace nikitina_v_hoar_sort_batcher {

namespace {

void QuickSortHoare(std::vector<int> &arr, int low, int high) {
  if (low >= high) {
    return;
  }
  std::vector<std::pair<int, int>> stack;
  stack.emplace_back(low, high);
  while (!stack.empty()) {
    auto [left_bound, right_bound] = stack.back();
    stack.pop_back();
    if (left_bound >= right_bound) {
      continue;
    }
    int pivot = arr[left_bound + ((right_bound - left_bound) / 2)];
    int left_idx = left_bound - 1;
    int right_idx = right_bound + 1;
    while (true) {
      left_idx++;
      while (arr[left_idx] < pivot) {
        left_idx++;
      }
      right_idx--;
      while (arr[right_idx] > pivot) {
        right_idx--;
      }
      if (left_idx >= right_idx) {
        break;
      }
      std::swap(arr[left_idx], arr[right_idx]);
    }
    stack.emplace_back(left_bound, right_idx);
    stack.emplace_back(right_idx + 1, right_bound);
  }
}

void ThreadedSort(std::vector<int> &arr, int num_threads) {
  int size = static_cast<int>(arr.size());
  if (size <= 1) {
    return;
  }
  int active_threads = std::min(num_threads, size / 2);
  if (active_threads <= 1) {
    QuickSortHoare(arr, 0, size - 1);
    return;
  }
  std::vector<std::thread> threads;
  int chunk_size = size / active_threads;
  for (int iter = 0; iter < active_threads; ++iter) {
    int start = iter * chunk_size;
    int end = (iter == active_threads - 1) ? size - 1 : (start + chunk_size - 1);
    threads.emplace_back([&arr, start, end]() { QuickSortHoare(arr, start, end); });
  }
  for (auto &thr : threads) {
    thr.join();
  }

  QuickSortHoare(arr, 0, size - 1);
}

void MpiCompareSwap(std::vector<int> &local_arr, int neighbor, bool keep_low) {
  int size = static_cast<int>(local_arr.size());
  std::vector<int> neighbor_arr(size);
  MPI_Sendrecv(local_arr.data(), size, MPI_INT, neighbor, 0, neighbor_arr.data(), size, MPI_INT, neighbor, 0,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  std::vector<int> full_merge(static_cast<size_t>(size) * 2);
  std::ranges::merge(local_arr, neighbor_arr, full_merge.begin());

  if (keep_low) {
    std::copy(full_merge.begin(), full_merge.begin() + size, local_arr.begin());
  } else {
    std::copy(full_merge.begin() + size, full_merge.end(), local_arr.begin());
  }
}

void BatcherExchange(std::vector<int> &local_data, int mpi_rank, int mpi_size, int p_step, int k_step) {
  for (int j_idx = k_step % p_step; j_idx + k_step < mpi_size; j_idx += (k_step << 1)) {
    for (int i_idx = 0; i_idx < std::min(k_step, mpi_size - j_idx - k_step); ++i_idx) {
      int r1 = j_idx + i_idx;
      int r2 = j_idx + i_idx + k_step;
      if (r1 / (p_step << 1) == r2 / (p_step << 1)) {
        if (mpi_rank == r1) {
          MpiCompareSwap(local_data, r2, true);
        } else if (mpi_rank == r2) {
          MpiCompareSwap(local_data, r1, false);
        }
      }
    }
  }
}

void BatcherMergeNetwork(std::vector<int> &local_data, int mpi_rank, int mpi_size) {
  for (int p_step = 1; p_step < mpi_size; p_step <<= 1) {
    for (int k_step = p_step; k_step > 0; k_step >>= 1) {
      BatcherExchange(local_data, mpi_rank, mpi_size, p_step, k_step);
    }
  }
}

}  // namespace

HoareSortBatcherALL::HoareSortBatcherALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool HoareSortBatcherALL::ValidationImpl() {
  return true;
}

bool HoareSortBatcherALL::PreProcessingImpl() {
  data_ = GetInput();
  return true;
}

bool HoareSortBatcherALL::RunImpl() {
  int mpi_rank = 0;
  int mpi_size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  int total_n = (mpi_rank == 0) ? static_cast<int>(data_.size()) : 0;
  MPI_Bcast(&total_n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (total_n == 0) {
    return true;
  }

  int chunk = (total_n + mpi_size - 1) / mpi_size;
  std::vector<int> local_data(chunk, std::numeric_limits<int>::max());
  std::vector<int> send_buffer;

  if (mpi_rank == 0) {
    send_buffer = data_;
    send_buffer.resize(static_cast<size_t>(chunk) * mpi_size, std::numeric_limits<int>::max());
  }

  MPI_Scatter(send_buffer.data(), chunk, MPI_INT, local_data.data(), chunk, MPI_INT, 0, MPI_COMM_WORLD);

  int hw_threads = static_cast<int>(std::thread::hardware_concurrency());
  ThreadedSort(local_data, hw_threads);

  BatcherMergeNetwork(local_data, mpi_rank, mpi_size);

  std::vector<int> gather_buffer;
  if (mpi_rank == 0) {
    gather_buffer.resize(static_cast<size_t>(chunk) * mpi_size);
  }

  MPI_Gather(local_data.data(), chunk, MPI_INT, gather_buffer.data(), chunk, MPI_INT, 0, MPI_COMM_WORLD);

  if (mpi_rank == 0) {
    gather_buffer.resize(total_n);
    data_ = std::move(gather_buffer);
  }

  return true;
}

bool HoareSortBatcherALL::PostProcessingImpl() {
  int mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  if (mpi_rank != 0) {
    data_.clear();
  }

  GetOutput() = data_;
  return true;
}

}  // namespace nikitina_v_hoar_sort_batcher
