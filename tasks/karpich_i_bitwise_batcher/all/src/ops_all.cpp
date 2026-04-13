#include "karpich_i_bitwise_batcher/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <random>
#include <utility>
#include <vector>

#include "karpich_i_bitwise_batcher/common/include/common.hpp"

namespace karpich_i_bitwise_batcher {

namespace {

int FindMaxParallel(const std::vector<int> &arr, int n) {
  int max_val = arr[0];
#pragma omp parallel for default(none) shared(arr, n) reduction(max : max_val)
  for (int i = 1; i < n; i++) {
    max_val = std::max(max_val, arr[i]);
  }
  return max_val;
}

void CountingPass(std::vector<int> &arr, std::vector<int> &buffer, int n, int shift) {
  std::vector<int> count(256, 0);
  int num_threads = omp_get_max_threads();
  std::vector<std::vector<int>> local_counts(num_threads, std::vector<int>(256, 0));

#pragma omp parallel default(none) shared(arr, shift, local_counts, n)
  {
    int tid = omp_get_thread_num();
#pragma omp for
    for (int i = 0; i < n; i++) {
      local_counts[tid][(arr[i] >> shift) & 0xFF]++;
    }
  }

  for (int ti = 0; ti < num_threads; ti++) {
    for (int i = 0; i < 256; i++) {
      count[i] += local_counts[ti][i];
    }
  }

  for (int i = 1; i < 256; i++) {
    count[i] += count[i - 1];
  }

  for (int i = n - 1; i >= 0; i--) {
    buffer[--count[(arr[i] >> shift) & 0xFF]] = arr[i];
  }
  arr = buffer;
}

void RadixSortPositive(std::vector<int> &arr) {
  int n = static_cast<int>(arr.size());
  if (n <= 1) {
    return;
  }

  int max_val = FindMaxParallel(arr, n);
  if (max_val == 0) {
    return;
  }

  std::vector<int> buffer(n);
  for (int shift = 0; shift < 32 && (max_val >> shift) > 0; shift += 8) {
    CountingPass(arr, buffer, n, shift);
  }
}

void RadixSort(std::vector<int> &arr) {
  int n = static_cast<int>(arr.size());
  if (n <= 1) {
    return;
  }

  std::vector<int> negative;
  std::vector<int> positive;
  for (int i = 0; i < n; i++) {
    if (arr[i] < 0) {
      negative.push_back(-arr[i]);
    } else {
      positive.push_back(arr[i]);
    }
  }

  RadixSortPositive(positive);
  RadixSortPositive(negative);

  int idx = 0;
  for (int i = static_cast<int>(negative.size()) - 1; i >= 0; i--) {
    arr[idx++] = -negative[i];
  }
  for (int x : positive) {
    arr[idx++] = x;
  }
}

struct MergeTask {
  int lo;
  int hi;
  int r;
};

std::vector<std::vector<std::pair<int, int>>> BuildMergeNetwork(int lo, int hi) {
  std::vector<std::vector<std::pair<int, int>>> levels;
  std::vector<MergeTask> current = {{lo, hi, 1}};

  while (!current.empty()) {
    std::vector<MergeTask> next;
    std::vector<std::pair<int, int>> comps;

    for (const auto &[tlo, thi, tr] : current) {
      int step = tr * 2;
      if (step < thi - tlo) {
        next.push_back({tlo, thi, step});
        next.push_back({tlo + tr, thi, step});
        for (int i = tlo + tr; i + tr <= thi; i += step) {
          comps.emplace_back(i, i + tr);
        }
      } else if (tlo + tr <= thi) {
        comps.emplace_back(tlo, tlo + tr);
      }
    }

    levels.push_back(std::move(comps));
    current = std::move(next);
  }

  return levels;
}

void ApplyComparatorNetwork(std::vector<int> &arr, const std::vector<std::vector<std::pair<int, int>>> &levels) {
  for (int lvl = static_cast<int>(levels.size()) - 1; lvl >= 0; lvl--) {
    const auto &level = levels[lvl];
    int level_size = static_cast<int>(level.size());
#pragma omp parallel for default(none) shared(arr, level, level_size)
    for (int i = 0; i < level_size; ++i) {
      int aa = level[i].first;
      int bb = level[i].second;
      if (arr[aa] > arr[bb]) {
        std::swap(arr[aa], arr[bb]);
      }
    }
  }
}

void BatcherMerge(std::vector<int> &arr, int lo, int hi) {
  auto levels = BuildMergeNetwork(lo, hi);
  ApplyComparatorNetwork(arr, levels);
}

void SortSingleProcess(std::vector<int> &data, int padded, int n) {
  int half = padded / 2;
  std::vector<int> left(data.begin(), data.begin() + half);
  std::vector<int> right(data.begin() + half, data.end());

  RadixSort(left);
  RadixSort(right);

  std::ranges::copy(left, data.begin());
  std::ranges::copy(right, data.begin() + half);

  BatcherMerge(data, 0, padded - 1);
  data.resize(n);
}

int PadToPowerOfTwo(std::vector<int> &data, int n) {
  if (n <= 1) {
    return n;
  }
  int padded = 1;
  while (padded < n) {
    padded *= 2;
  }
  int max_elem = *std::ranges::max_element(data);
  data.resize(padded, max_elem);
  return padded;
}

int ComputeNumTasks(int num_ranks, int padded) {
  int num_tasks = 1;
  while (num_tasks * 2 <= num_ranks && num_tasks * 2 <= padded) {
    num_tasks *= 2;
  }
  return num_tasks;
}

void MergeChunks(std::vector<int> &arr, int chunk_size, int padded) {
  for (int step = chunk_size * 2; step <= padded; step *= 2) {
#pragma omp parallel for default(none) shared(step, padded, arr)
    for (int i = 0; i < padded; i += step) {
      BatcherMerge(arr, i, i + step - 1);
    }
  }
}

}  // namespace

KarpichIBitwiseBatcherALL::KarpichIBitwiseBatcherALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool KarpichIBitwiseBatcherALL::ValidationImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    return GetInput() > 0;
  }
  return true;
}

bool KarpichIBitwiseBatcherALL::PreProcessingImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    int n = GetInput();
    data_.resize(n);
    std::mt19937 gen(static_cast<unsigned int>(n));
    std::uniform_int_distribution<int> dist(-1000, 1000);
    for (int i = 0; i < n; i++) {
      data_[i] = dist(gen);
    }
  }
  return true;
}

bool KarpichIBitwiseBatcherALL::RunImpl() {
  int rank = 0;
  int num_ranks = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

  int n = 0;
  int padded = 1;

  if (rank == 0) {
    n = static_cast<int>(data_.size());
    padded = PadToPowerOfTwo(data_, n);
  }

  MPI_Bcast(&padded, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (padded <= 1) {
    return true;
  }

  int num_tasks = ComputeNumTasks(num_ranks, padded);

  if (num_tasks == 1) {
    if (rank == 0) {
      SortSingleProcess(data_, padded, n);
    }
    return true;
  }

  MPI_Comm active_comm = MPI_COMM_NULL;
  int color = (rank < num_tasks) ? 1 : MPI_UNDEFINED;
  MPI_Comm_split(MPI_COMM_WORLD, color, rank, &active_comm);

  if (active_comm != MPI_COMM_NULL) {
    int chunk_size = padded / num_tasks;
    std::vector<int> local_data(chunk_size);

    MPI_Scatter(data_.data(), chunk_size, MPI_INT, local_data.data(), chunk_size, MPI_INT, 0, active_comm);

    RadixSort(local_data);

    MPI_Gather(local_data.data(), chunk_size, MPI_INT, data_.data(), chunk_size, MPI_INT, 0, active_comm);

    if (rank == 0) {
      MergeChunks(data_, chunk_size, padded);
      data_.resize(n);
    }
    MPI_Comm_free(&active_comm);
  }

  return true;
}

bool KarpichIBitwiseBatcherALL::PostProcessingImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    for (int i = 1; std::cmp_less(i, data_.size()); i++) {
      if (data_[i] < data_[i - 1]) {
        return false;
      }
    }
  }
  GetOutput() = GetInput();
  return true;
}

}  // namespace karpich_i_bitwise_batcher
