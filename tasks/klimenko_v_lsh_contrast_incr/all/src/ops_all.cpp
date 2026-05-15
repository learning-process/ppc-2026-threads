#include "klimenko_v_lsh_contrast_incr/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <climits>
#include <cstddef>
#include <thread>
#include <utility>
#include <vector>

#include "klimenko_v_lsh_contrast_incr/common/include/common.hpp"
#include "oneapi/tbb/parallel_for.h"
#include "util/include/util.hpp"

namespace klimenko_v_lsh_contrast_incr {

namespace {

void GetThreadRange(size_t tid, size_t total, size_t num_t, size_t &begin, size_t &end) {
  size_t chunk = total / num_t;
  begin = tid * chunk;
  end = (tid == num_t - 1) ? total : begin + chunk;
}

std::pair<int, int> FindMinMaxSTL(const std::vector<int> &data, int num_threads) {
  const size_t size = data.size();
  if (size == 0) {
    return {INT_MAX, INT_MIN};
  }

  std::vector<std::thread> threads(num_threads);
  std::vector<int> local_min(num_threads);
  std::vector<int> local_max(num_threads);

  for (int tid = 0; tid < num_threads; tid++) {
    threads[tid] = std::thread([&, tid]() {
      size_t begin = 0;
      size_t end = 0;
      GetThreadRange(tid, size, num_threads, begin, end);

      if (begin >= size) {
        local_min[tid] = INT_MAX;
        local_max[tid] = INT_MIN;
        return;
      }

      int min_val = data[begin];
      int max_val = data[begin];

      for (size_t i = begin; i < end; i++) {
        min_val = std::min(min_val, data[i]);
        max_val = std::max(max_val, data[i]);
      }

      local_min[tid] = min_val;
      local_max[tid] = max_val;
    });
  }

  for (auto &th : threads) {
    th.join();
  }

  int global_min = *std::ranges::min_element(local_min);
  int global_max = *std::ranges::max_element(local_max);

  return {global_min, global_max};
}

void StretchContrastSTL(const std::vector<int> &input, std::vector<int> &output, int min_val, int max_val,
                        int num_threads) {
  const size_t size = input.size();
  if (size == 0) {
    return;
  }

  std::vector<std::thread> threads(num_threads);

  for (int tid = 0; tid < num_threads; tid++) {
    threads[tid] = std::thread([&, tid]() {
      size_t begin = 0;
      size_t end = 0;
      GetThreadRange(tid, size, num_threads, begin, end);

      for (size_t i = begin; i < end; i++) {
        output[i] = ((input[i] - min_val) * 255) / (max_val - min_val);
      }
    });
  }

  for (auto &th : threads) {
    th.join();
  }
}

}  // namespace

KlimenkoVLSHContrastIncrALL::KlimenkoVLSHContrastIncrALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool KlimenkoVLSHContrastIncrALL::ValidationImpl() {
  int initialized = 0;
  MPI_Initialized(&initialized);
  if (initialized != 0) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
      return !GetInput().empty();
    }
    return true;
  }
  return !GetInput().empty();
}

bool KlimenkoVLSHContrastIncrALL::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return true;
}

bool KlimenkoVLSHContrastIncrALL::RunImpl() {
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const auto &input = GetInput();
  auto &output = GetOutput();

  int total_size = 0;
  if (rank == 0) {
    total_size = static_cast<int>(input.size());
  }
  MPI_Bcast(&total_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (total_size == 0) {
    return false;
  }

  const int num_stl_threads = ppc::util::GetNumThreads();

  std::vector<int> recv_counts(size);
  std::vector<int> displs(size);

  int chunk = total_size / size;
  int remainder = total_size % size;
  for (int i = 0; i < size; i++) {
    recv_counts[i] = chunk + (i < remainder ? 1 : 0);
    displs[i] = (i == 0) ? 0 : displs[i - 1] + recv_counts[i - 1];
  }

  int local_size = recv_counts[rank];
  std::vector<int> local_data(local_size);

  MPI_Scatterv(input.data(), recv_counts.data(), displs.data(), MPI_INT, local_data.data(), local_size, MPI_INT, 0,
               MPI_COMM_WORLD);

  auto [local_min, local_max] = FindMinMaxSTL(local_data, num_stl_threads);

  int global_min = 0;
  int global_max = 0;
  MPI_Allreduce(&local_min, &global_min, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&local_max, &global_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  std::vector<int> local_output(local_size);

  if (global_max == global_min) {
    std::ranges::copy(local_data, local_output.begin());
  } else {
    StretchContrastSTL(local_data, local_output, global_min, global_max, num_stl_threads);
  }

  std::vector<int> recv_buffer;
  if (rank == 0) {
    recv_buffer.resize(total_size);
  }

  MPI_Gatherv(local_output.data(), local_size, MPI_INT, recv_buffer.data(), recv_counts.data(), displs.data(), MPI_INT,
              0, MPI_COMM_WORLD);

  if (rank == 0) {
    std::ranges::copy(recv_buffer, output.begin());
  }

  return true;
}

bool KlimenkoVLSHContrastIncrALL::PostProcessingImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    return GetOutput().size() == GetInput().size();
  }
  return true;
}

}  // namespace klimenko_v_lsh_contrast_incr
