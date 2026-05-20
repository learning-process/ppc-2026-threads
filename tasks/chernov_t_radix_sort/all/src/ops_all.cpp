#include "chernov_t_radix_sort/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "chernov_t_radix_sort/common/include/common.hpp"

namespace chernov_t_radix_sort {

ChernovTRadixSortALL::ChernovTRadixSortALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool ChernovTRadixSortALL::ValidationImpl() {
  return true;
}

bool ChernovTRadixSortALL::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

constexpr int kBitsPerDigit = 8;
constexpr int kRadix = 1 << kBitsPerDigit;
constexpr uint32_t kSignMask = 0x80000000U;

void ChernovTRadixSortALL::RadixSortLSDParallelOMP(std::vector<int> &data) {
  if (data.empty()) {
    return;
  }
  const size_t n = data.size();

  std::vector<uint32_t> temp(n);
#pragma omp parallel for schedule(static) default(none) shared(data, temp, n)
  for (size_t i = 0; i < n; ++i) {
    temp[i] = static_cast<uint32_t>(data[i]) ^ kSignMask;
  }

  std::vector<uint32_t> buffer(n);
  int num_threads = omp_get_max_threads();

  for (int byte = 0; byte < 4; ++byte) {
    const int shift = byte * kBitsPerDigit;

    std::vector<std::vector<int>> local_counts(static_cast<size_t>(num_threads), std::vector<int>(kRadix, 0));

#pragma omp parallel for schedule(static) default(none) shared(temp, local_counts, n, shift)
    for (size_t i = 0; i < n; ++i) {
      int thread_idx = omp_get_thread_num();
      int digit = static_cast<int>((temp[i] >> shift) & 0xFFU);
      local_counts[static_cast<size_t>(thread_idx)][static_cast<size_t>(digit)]++;
    }

    std::vector<int> global_count(kRadix, 0);
    for (int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
      for (int digit_idx = 0; digit_idx < kRadix; ++digit_idx) {
        global_count[digit_idx] += local_counts[static_cast<size_t>(thread_idx)][static_cast<size_t>(digit_idx)];
      }
    }

    for (int i = 1; i < kRadix; ++i) {
      global_count[i] += global_count[i - 1];
    }

    std::vector<int> count = global_count;
    for (size_t i = n; i-- > 0;) {
      int digit = static_cast<int>((temp[i] >> shift) & 0xFFU);
      buffer[static_cast<size_t>(--count[static_cast<size_t>(digit)])] = temp[i];
    }

    temp.swap(buffer);
  }

#pragma omp parallel for schedule(static) default(none) shared(data, temp, n)
  for (size_t i = 0; i < n; ++i) {
    data[i] = static_cast<int>(temp[i] ^ kSignMask);
  }
}

void ChernovTRadixSortALL::SimpleMerge(const std::vector<int> &left, const std::vector<int> &right,
                                       std::vector<int> &result) {
  result.resize(left.size() + right.size());
  std::ranges::merge(left, right, result.begin());
}

void ChernovTRadixSortALL::ComputeChunkSizes(int num_processes, size_t total_elements, std::vector<int> &recv_counts,
                                             std::vector<int> &displs) {
  int base = static_cast<int>(total_elements / static_cast<size_t>(num_processes));
  int remainder = static_cast<int>(total_elements % static_cast<size_t>(num_processes));
  int current_disp = 0;
  for (int i = 0; i < num_processes; ++i) {
    recv_counts[i] = base + (i < remainder ? 1 : 0);
    displs[i] = current_disp;
    current_disp += recv_counts[i];
  }
}

void ChernovTRadixSortALL::MergeChunksOnRank0(const std::vector<int> &global_result,
                                              const std::vector<int> &recv_counts, const std::vector<int> &displs,
                                              std::vector<int> &output) {
  if (global_result.empty()) {
    output.clear();
    return;
  }

  int num_processes = static_cast<int>(recv_counts.size());

  int start_idx = -1;
  int offset = 0;
  for (int i = 0; i < num_processes; ++i) {
    if (recv_counts[i] > 0) {
      start_idx = i;
      offset = displs[i];
      break;
    }
  }

  if (start_idx == -1) {
    output.clear();
    return;
  }

  std::vector<int> merged(global_result.begin() + offset, global_result.begin() + offset + recv_counts[start_idx]);
  offset += recv_counts[start_idx];

  for (int i = start_idx + 1; i < num_processes; ++i) {
    if (recv_counts[i] > 0) {
      std::vector<int> next_part(global_result.begin() + offset, global_result.begin() + offset + recv_counts[i]);
      std::vector<int> new_merged;
      SimpleMerge(merged, next_part, new_merged);
      merged = std::move(new_merged);
      offset += recv_counts[i];
    }
  }

  output = std::move(merged);
}

bool ChernovTRadixSortALL::RunImpl() {
  auto &input_data = GetInput();
  int current_rank = -1;
  int num_processes = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &current_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_processes);

  const size_t total_elements = input_data.size();

  if (total_elements == 0) {
    GetOutput().clear();
    return true;
  }

  std::vector<int> recv_counts(num_processes, 0);
  std::vector<int> displs(num_processes, 0);

  if (current_rank == 0) {
    ComputeChunkSizes(num_processes, total_elements, recv_counts, displs);
  }

  int local_n = 0;
  MPI_Scatter(recv_counts.data(), 1, MPI_INT, &local_n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  std::vector<int> local_data;
  if (local_n > 0) {
    local_data.resize(static_cast<size_t>(local_n));
    MPI_Scatterv(input_data.data(), recv_counts.data(), displs.data(), MPI_INT, local_data.data(), local_n, MPI_INT, 0,
                 MPI_COMM_WORLD);
    RadixSortLSDParallelOMP(local_data);
  }

  std::vector<int> global_result;
  if (current_rank == 0) {
    global_result.resize(total_elements);
  }

  MPI_Gatherv(local_data.empty() ? nullptr : local_data.data(), local_n, MPI_INT, global_result.data(),
              recv_counts.data(), displs.data(), MPI_INT, 0, MPI_COMM_WORLD);

  if (current_rank == 0) {
    MergeChunksOnRank0(global_result, recv_counts, displs, GetOutput());
  } else {
    GetOutput().clear();
  }

  int out_size = static_cast<int>(GetOutput().size());
  MPI_Bcast(&out_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (current_rank != 0) {
    GetOutput().resize(static_cast<size_t>(out_size));
  }

  if (out_size > 0) {
    MPI_Bcast(GetOutput().data(), out_size, MPI_INT, 0, MPI_COMM_WORLD);
  }

  return true;
}

bool ChernovTRadixSortALL::PostProcessingImpl() {
  return std::is_sorted(GetOutput().begin(), GetOutput().end());
}

}  // namespace chernov_t_radix_sort
