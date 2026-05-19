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

void MergeBlocksLocal(std::vector<int> &arr, int start1, int start2, int chunk_size) {
  std::vector<int> buffer(static_cast<std::size_t>(chunk_size) * 2);
  int i = start1;
  int j = start2;
  int k = 0;

  int end1 = start1 + chunk_size;
  int end2 = start2 + chunk_size;

  while (i < end1 && j < end2) {
    if (arr[i] <= arr[j]) {
      buffer[k++] = arr[i++];
    } else {
      buffer[k++] = arr[j++];
    }
  }

  while (i < end1) {
    buffer[k++] = arr[i++];
  }
  while (j < end2) {
    buffer[k++] = arr[j++];
  }

  for (int idx = 0; idx < chunk_size; ++idx) {
    arr[start1 + idx] = buffer[idx];
    arr[start2 + idx] = buffer[chunk_size + idx];
  }
}

void ProcessMergeStepLocal(std::vector<int> &data, int step_j, int step_p, int step_k, int chunk_size) {
  for (int i = 0; i < step_k; ++i) {
    if ((step_j + i) / (step_p * 2) == (step_j + i + step_k) / (step_p * 2)) {
      int start_a = (step_j + i) * chunk_size;
      int start_b = (step_j + i + step_k) * chunk_size;
      MergeBlocksLocal(data, start_a, start_b, chunk_size);
    }
  }
}

void BatcherMergeLocalOMP(std::vector<int> &data, int num_threads, int chunk_size) {
  for (int step_p = 1; step_p < num_threads; step_p *= 2) {
    for (int step_k = step_p; step_k > 0; step_k /= 2) {
      int limit = num_threads - step_k;
      int start_j = step_k % step_p;

      if (start_j <= limit) {
#pragma omp parallel for default(none) shared(data, num_threads, chunk_size, step_p, step_k, start_j, limit)
        for (int step_j = start_j; step_j <= limit; step_j += (step_k * 2)) {
          ProcessMergeStepLocal(data, step_j, step_p, step_k, chunk_size);
        }
      }
    }
  }
}

void FindPartnerForMPI(int rank, int active_size, int step_p, int step_k, int &partner, bool &is_left) {
  partner = -1;
  is_left = false;
  for (int step_j = step_k % step_p; step_j + step_k < active_size; step_j += (step_k * 2)) {
    for (int i = 0; i < std::min(step_k, active_size - step_j - step_k); ++i) {
      if ((step_j + i) / (step_p * 2) == (step_j + i + step_k) / (step_p * 2)) {
        int r1 = step_j + i;
        int r2 = step_j + i + step_k;
        if (rank == r1) {
          partner = r2;
          is_left = true;
          return;
        }
        if (rank == r2) {
          partner = r1;
          is_left = false;
          return;
        }
      }
    }
  }
}

void BatcherMergeGlobalMPI(std::vector<int> &local_data, int rank, int active_size) {
  int chunk_size = static_cast<int>(local_data.size());

  for (int step_p = 1; step_p < active_size; step_p *= 2) {
    for (int step_k = step_p; step_k > 0; step_k /= 2) {
      int partner = -1;
      bool is_left = false;

      FindPartnerForMPI(rank, active_size, step_p, step_k, partner, is_left);

      if (partner != -1) {
        std::vector<int> recv_data(static_cast<std::size_t>(chunk_size));
        MPI_Status status;
        MPI_Sendrecv(local_data.data(), chunk_size, MPI_INT, partner, 0, recv_data.data(), chunk_size, MPI_INT, partner,
                     0, MPI_COMM_WORLD, &status);

        std::vector<int> merged(static_cast<std::size_t>(chunk_size) * 2);
        std::ranges::merge(local_data, recv_data, merged.begin());

        if (is_left) {
          std::copy(merged.begin(), merged.begin() + chunk_size, local_data.begin());
        } else {
          std::copy(merged.begin() + chunk_size, merged.end(), local_data.begin());
        }
      }
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
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    input_ = GetInput();
  }
  return true;
}

bool ShekhirevHoareBatcherSortALL::RunImpl() {
  int rank = 0;
  int world_size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  int orig_size = 0;
  if (rank == 0) {
    orig_size = static_cast<int>(input_.size());
  }
  MPI_Bcast(&orig_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (orig_size <= 1) {
    if (rank == 0) {
      output_ = input_;
    }
    return true;
  }

  int active_size = 1;
  while (active_size * 2 <= world_size && active_size * 2 <= orig_size) {
    active_size *= 2;
  }

  int mpi_padding = (active_size - (orig_size % active_size)) % active_size;
  int mpi_chunk_size = (orig_size + mpi_padding) / active_size;

  std::vector<int> counts(static_cast<std::size_t>(world_size), 0);
  std::vector<int> displs(static_cast<std::size_t>(world_size), 0);
  for (int i = 0; i < active_size; ++i) {
    counts[static_cast<std::size_t>(i)] = mpi_chunk_size;
    displs[static_cast<std::size_t>(i)] = i * mpi_chunk_size;
  }

  std::vector<int> global_data;
  if (rank == 0) {
    global_data = input_;
    if (mpi_padding > 0) {
      global_data.insert(global_data.end(), mpi_padding, std::numeric_limits<int>::max());
    }
  }

  std::vector<int> local_data(rank < active_size ? static_cast<std::size_t>(mpi_chunk_size) : 0);
  const int *send_buf = (rank == 0) ? global_data.data() : nullptr;

  MPI_Scatterv(send_buf, counts.data(), displs.data(), MPI_INT, local_data.data(),
               counts[static_cast<std::size_t>(rank)], MPI_INT, 0, MPI_COMM_WORLD);

  if (rank < active_size && mpi_chunk_size > 0) {
    int max_threads = omp_get_max_threads();
    if (max_threads <= 0) {
      max_threads = 1;
    }

    int num_omp_threads = 1;
    while (num_omp_threads * 2 <= max_threads && num_omp_threads * 2 <= mpi_chunk_size) {
      num_omp_threads *= 2;
    }

    if (num_omp_threads == 1) {
      OptimizedHoareSort(local_data, 0, mpi_chunk_size - 1);
    } else {
      int local_pad = (num_omp_threads - (mpi_chunk_size % num_omp_threads)) % num_omp_threads;
      if (local_pad > 0) {
        local_data.insert(local_data.end(), local_pad, std::numeric_limits<int>::max());
      }
      int local_chunk = static_cast<int>(local_data.size()) / num_omp_threads;

#pragma omp parallel for default(none) shared(local_data, num_omp_threads, local_chunk)
      for (int i = 0; i < num_omp_threads; ++i) {
        OptimizedHoareSort(local_data, i * local_chunk, ((i + 1) * local_chunk) - 1);
      }
      BatcherMergeLocalOMP(local_data, num_omp_threads, local_chunk);

      if (local_pad > 0) {
        local_data.resize(static_cast<std::size_t>(mpi_chunk_size));
      }
    }

    BatcherMergeGlobalMPI(local_data, rank, active_size);
  }

  int *recv_buf = (rank == 0) ? global_data.data() : nullptr;
  MPI_Gatherv(local_data.data(), counts[static_cast<std::size_t>(rank)], MPI_INT, recv_buf, counts.data(),
              displs.data(), MPI_INT, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    output_ = std::move(global_data);
  }

  return true;
}

bool ShekhirevHoareBatcherSortALL::PostProcessingImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    output_.resize(input_.size());
    GetOutput() = std::move(output_);
  }
  return true;
}

}  // namespace shekhirev_v_hoare_batcher_sort
