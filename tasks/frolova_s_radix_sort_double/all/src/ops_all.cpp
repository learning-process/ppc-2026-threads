#include "frolova_s_radix_sort_double/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>
#include <oneapi/tbb/parallel_for.h>

#include <algorithm>
#include <atomic>
#include <bit>
#include <cmath>
#include <cstdint>
#include <thread>
#include <vector>

#include "frolova_s_radix_sort_double/common/include/common.hpp"

namespace frolova_s_radix_sort_double {

std::vector<double> FrolovaSRadixSortDoubleALL::MergeSorted(const std::vector<double> &a,
                                                            const std::vector<double> &b) {
  std::vector<double> result;
  result.reserve(a.size() + b.size());
  std::merge(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(result));
  return result;
}

FrolovaSRadixSortDoubleALL::FrolovaSRadixSortDoubleALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool FrolovaSRadixSortDoubleALL::ValidationImpl() {
  return !GetInput().empty();
}

bool FrolovaSRadixSortDoubleALL::PreProcessingImpl() {
  return true;
}

bool FrolovaSRadixSortDoubleALL::RunImpl() {
  int rank;
  int size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  std::vector<double> global_input;
  if (rank == 0) {
    global_input = GetInput();
  }

  int total_n = 0;
  if (rank == 0) {
    total_n = static_cast<int>(global_input.size());
  }
  MPI_Bcast(&total_n, 1, MPI_INT, 0, MPI_COMM_WORLD);

  int chunk_size = total_n / size;
  int remainder = total_n % size;

  std::vector<int> send_counts(size, chunk_size);
  for (int i = 0; i < remainder; ++i) {
    send_counts[i]++;
  }

  std::vector<int> displs(size, 0);
  for (int i = 1; i < size; ++i) {
    displs[i] = displs[i - 1] + send_counts[i - 1];
  }

  std::vector<double> local_data(send_counts[rank]);
  MPI_Scatterv(global_input.data(), send_counts.data(), displs.data(), MPI_DOUBLE, local_data.data(), send_counts[rank],
               MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (!local_data.empty()) {
    const int kRadix = 256;
    const int kNumPasses = 8;
    std::vector<double> temp(local_data.size());

    for (int pass = 0; pass < kNumPasses; ++pass) {
      std::vector<std::atomic<int>> counts(kRadix);
      for (int i = 0; i < kRadix; ++i) {
        counts[i].store(0);
      }

      tbb::parallel_for(tbb::blocked_range<size_t>(0, local_data.size()), [&](const tbb::blocked_range<size_t> &r) {
        for (size_t i = r.begin(); i < r.end(); ++i) {
          uint64_t bits = std::bit_cast<uint64_t>(local_data[i]);
          int byte = (bits >> (pass * 8)) & 0xFF;
          counts[byte].fetch_add(1, std::memory_order_relaxed);
        }
      });

      std::vector<int> offsets(kRadix);
      int total = 0;
      for (int i = 0; i < kRadix; ++i) {
        offsets[i] = total;
        total += counts[i].load();
      }

      std::atomic<int> sync_counter(0);
      std::thread sync_thread([&sync_counter]() { sync_counter.fetch_add(1); });
      sync_thread.join();

      std::vector<int> current_offsets = offsets;
      for (double val : local_data) {
        uint64_t bits = std::bit_cast<uint64_t>(val);
        int byte = (bits >> (pass * 8)) & 0xFF;
        temp[current_offsets[byte]++] = val;
      }
      local_data.swap(temp);
    }

    std::vector<double> negatives;
    std::vector<double> positives;

#pragma omp parallel
    {
      std::vector<double> neg_local, pos_local;
#pragma omp for nowait
      for (int i = 0; i < static_cast<int>(local_data.size()); ++i) {
        if (local_data[i] < 0) {
          neg_local.push_back(local_data[i]);
        } else {
          pos_local.push_back(local_data[i]);
        }
      }
#pragma omp critical
      {
        negatives.insert(negatives.end(), neg_local.begin(), neg_local.end());
        positives.insert(positives.end(), pos_local.begin(), pos_local.end());
      }
    }
    std::reverse(negatives.begin(), negatives.end());

    local_data.clear();
    local_data.insert(local_data.end(), negatives.begin(), negatives.end());
    local_data.insert(local_data.end(), positives.begin(), positives.end());
  }

  if (rank == 0) {
    std::vector<std::vector<double>> all_parts(size);
    all_parts[0] = local_data;
    for (int i = 1; i < size; ++i) {
      all_parts[i].resize(send_counts[i]);
      MPI_Recv(all_parts[i].data(), send_counts[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    std::vector<double> merged_result = all_parts[0];
    for (int i = 1; i < size; ++i) {
      merged_result = MergeSorted(merged_result, all_parts[i]);
    }
    GetOutput() = std::move(merged_result);
  } else {
    MPI_Send(local_data.data(), static_cast<int>(local_data.size()), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }

  return true;
}

bool FrolovaSRadixSortDoubleALL::PostProcessingImpl() {
  return true;
}

}  // namespace frolova_s_radix_sort_double
