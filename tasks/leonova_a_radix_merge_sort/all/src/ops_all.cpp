#include "leonova_a_radix_merge_sort/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <future>
#include <utility>
#include <vector>

#include "leonova_a_radix_merge_sort/common/include/common.hpp"
#include "util/include/util.hpp"

static_assert(sizeof(int64_t) == sizeof(std::int64_t));

namespace leonova_a_radix_merge_sort {

LeonovaARadixMergeSortALL::LeonovaARadixMergeSortALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<int64_t>(GetInput().size());
}

bool LeonovaARadixMergeSortALL::ValidationImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int is_valid = 1;
  if (rank == 0) {
    is_valid = GetInput().empty() ? 0 : 1;
  }

  MPI_Bcast(&is_valid, 1, MPI_INT, 0, MPI_COMM_WORLD);
  return is_valid == 1;
}

bool LeonovaARadixMergeSortALL::PreProcessingImpl() {
  return true;
}

bool LeonovaARadixMergeSortALL::PostProcessingImpl() {
  return true;
}

namespace {

void MpiSendInt64(const std::vector<int64_t> &buf, int count, int dest, int tag) {
  MPI_Send(static_cast<const void *>(buf.data()), count, MPI_INT64_T, dest, tag, MPI_COMM_WORLD);
}

void MpiRecvInt64(std::vector<int64_t> &buf, int count, int src, int tag) {
  MPI_Recv(static_cast<void *>(buf.data()), count, MPI_INT64_T, src, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void MpiBcastInt64(std::vector<int64_t> &buf, int count) {
  MPI_Bcast(static_cast<void *>(buf.data()), count, MPI_INT64_T, 0, MPI_COMM_WORLD);
}

void MpiScattervInt64(const std::vector<int64_t> &sendbuf, const std::vector<int> &counts,
                      const std::vector<int> &displs, std::vector<int64_t> &recvbuf, int recv_count, int root) {
  const void *send_ptr = static_cast<const void *>(sendbuf.data());
  void *recv_ptr = static_cast<void *>(recvbuf.data());
  MPI_Scatterv(send_ptr, counts.data(), displs.data(), MPI_INT64_T, recv_ptr, recv_count, MPI_INT64_T, root,
               MPI_COMM_WORLD);
}

}  // namespace

inline uint64_t LeonovaARadixMergeSortALL::ToUnsignedValue(int64_t value) {
  return static_cast<uint64_t>(value) ^ kSignBitMask;
}

void LeonovaARadixMergeSortALL::FillKeys(const std::vector<int64_t> &arr, size_t left, size_t size,
                                         std::vector<uint64_t> &keys) {
  for (size_t i = 0; i < size; ++i) {
    keys[i] = ToUnsignedValue(arr[left + i]);
  }
}

void LeonovaARadixMergeSortALL::CountBytes(const std::vector<uint64_t> &keys, size_t begin, size_t end, int shift,
                                           std::vector<size_t> &counts) {
  std::ranges::fill(counts, 0);
  for (size_t i = begin; i < end; ++i) {
    ++counts[(keys[i] >> shift) & 0xFFU];
  }
}

void LeonovaARadixMergeSortALL::ScatterBytes(const std::vector<uint64_t> &src_keys, const std::vector<int64_t> &src_arr,
                                             size_t left, size_t begin, size_t end, int shift,
                                             std::vector<size_t> &offsets, std::vector<int64_t> &dst_arr,
                                             std::vector<uint64_t> &dst_keys) {
  for (size_t i = begin; i < end; ++i) {
    const size_t bucket = (src_keys[i] >> shift) & 0xFFU;
    const size_t pos = offsets[bucket]++;
    dst_arr[pos] = src_arr[left + i];
    dst_keys[pos] = src_keys[i];
  }
}

void LeonovaARadixMergeSortALL::ParallelRadixSort(std::vector<int64_t> &arr, size_t left, size_t right) {
  const size_t size = right - left;
  if (size <= 1) {
    return;
  }

  const auto hw = static_cast<size_t>(std::max(1, ppc::util::GetNumThreads()));
  const size_t num_threads = std::min(hw, size);
  const size_t chunk = (size + num_threads - 1) / num_threads;

  std::vector<size_t> boundaries(num_threads + 1);
  for (size_t tndex = 0; tndex <= num_threads; ++tndex) {
    boundaries[tndex] = left + std::min(tndex * chunk, size);
  }

  auto sort_chunk = [&](size_t cl, size_t cr) {
    const size_t sz = cr - cl;
    std::vector<uint64_t> keys(sz);
    std::vector<int64_t> tmp_arr(sz);
    std::vector<uint64_t> tmp_keys(sz);
    std::vector<size_t> counts(kNumCounters);
    std::vector<size_t> offsets(kNumCounters);

    FillKeys(arr, cl, sz, keys);

    for (int byte_pos = 0; byte_pos < kNumBytes; ++byte_pos) {
      const int shift = byte_pos * kByteSize;

      CountBytes(keys, 0, sz, shift, counts);

      size_t sum = 0;
      for (size_t bndex = 0; bndex < kNumCounters; ++bndex) {
        offsets[bndex] = sum;
        sum += counts[bndex];
      }

      ScatterBytes(keys, arr, cl, 0, sz, shift, offsets, tmp_arr, tmp_keys);

      std::ranges::copy(tmp_arr, arr.begin() + static_cast<std::ptrdiff_t>(cl));
      keys.swap(tmp_keys);
    }
  };

  std::vector<std::future<void>> futures;
  futures.reserve(num_threads);
  for (size_t tndex = 0; tndex < num_threads; ++tndex) {
    futures.push_back(std::async(std::launch::async, sort_chunk, boundaries[tndex], boundaries[tndex + 1]));
  }
  for (auto &f : futures) {
    f.get();
  }

  MergeChunks(arr, boundaries);
}

void LeonovaARadixMergeSortALL::SimpleMerge(std::vector<int64_t> &arr, size_t left, size_t mid, size_t right) {
  std::vector<int64_t> merged(right - left);

  size_t i = left;
  size_t j = mid;
  size_t k = 0;

  while (i < mid && j < right) {
    merged[k++] = (arr[i] <= arr[j]) ? arr[i++] : arr[j++];
  }
  while (i < mid) {
    merged[k++] = arr[i++];
  }
  while (j < right) {
    merged[k++] = arr[j++];
  }

  std::ranges::copy(merged, arr.begin() + static_cast<std::ptrdiff_t>(left));
}

void LeonovaARadixMergeSortALL::MergeChunks(std::vector<int64_t> &arr, const std::vector<size_t> &boundaries) {
  const size_t num_chunks = boundaries.size() - 1;
  if (num_chunks <= 1) {
    return;
  }

  size_t step = 1;
  while (step < num_chunks) {
    const size_t double_step = step * 2;
    size_t t = 0;
    while (t < num_chunks) {
      const size_t mid_idx = std::min(t + step, num_chunks);
      const size_t right_idx = std::min(t + double_step, num_chunks);
      if (mid_idx < right_idx) {
        SimpleMerge(arr, boundaries[t], boundaries[mid_idx], boundaries[right_idx]);
      }
      t += double_step;
    }
    step = double_step;
  }
}

void LeonovaARadixMergeSortALL::MergeSortedVectors(std::vector<int64_t> &local, const std::vector<int64_t> &incoming) {
  std::vector<int64_t> merged;
  merged.reserve(local.size() + incoming.size());

  size_t i = 0;
  size_t j = 0;

  while (i < local.size() && j < incoming.size()) {
    merged.push_back((local[i] <= incoming[j]) ? local[i++] : incoming[j++]);
  }
  while (i < local.size()) {
    merged.push_back(local[i++]);
  }
  while (j < incoming.size()) {
    merged.push_back(incoming[j++]);
  }

  local = std::move(merged);
}

void LeonovaARadixMergeSortALL::HierarchicalMerge(std::vector<int64_t> &local_data, int rank, int world_size) {
  int step = 1;

  while (step < world_size) {
    const int mask = step * 2;

    if ((rank & (mask - 1)) == 0) {
      const int partner = rank + step;
      if (partner < world_size) {
        int incoming_size = 0;
        MPI_Recv(&incoming_size, 1, MPI_INT, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if (incoming_size > 0) {
          std::vector<int64_t> incoming(static_cast<size_t>(incoming_size));
          MpiRecvInt64(incoming, incoming_size, partner, 1);
          MergeSortedVectors(local_data, incoming);
        }
      }
    } else if ((rank & (mask - 1)) == step) {
      const int partner = rank - step;
      const int send_size = static_cast<int>(local_data.size());
      MPI_Send(&send_size, 1, MPI_INT, partner, 0, MPI_COMM_WORLD);

      if (send_size > 0) {
        MpiSendInt64(local_data, send_size, partner, 1);
      }

      local_data.clear();
      local_data.shrink_to_fit();
    }

    MPI_Barrier(MPI_COMM_WORLD);
    step *= 2;
  }
}

void LeonovaARadixMergeSortALL::BroadcastResult(std::vector<int64_t> &local_data, int rank, int total) {
  int result_size = (rank == 0) ? static_cast<int>(local_data.size()) : 0;
  MPI_Bcast(&result_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    local_data.resize(static_cast<size_t>(total));
  }

  MpiBcastInt64(local_data, result_size);
}

bool LeonovaARadixMergeSortALL::RunImpl() {
  if (!ValidationImpl()) {
    return false;
  }

  int rank = 0;
  int world_size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  const auto total = static_cast<int>(GetInput().size());

  std::vector<int> send_counts(static_cast<size_t>(world_size), 0);
  std::vector<int> send_displs(static_cast<size_t>(world_size), 0);

  if (rank == 0) {
    const int base = total / world_size;
    const int remainder = total % world_size;
    for (int pndex = 0; pndex < world_size; ++pndex) {
      send_counts[static_cast<size_t>(pndex)] = base + (pndex < remainder ? 1 : 0);
    }
    int offset = 0;
    for (int pndex = 0; pndex < world_size; ++pndex) {
      send_displs[static_cast<size_t>(pndex)] = offset;
      offset += send_counts[static_cast<size_t>(pndex)];
    }
  }

  MPI_Bcast(send_counts.data(), world_size, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(send_displs.data(), world_size, MPI_INT, 0, MPI_COMM_WORLD);

  const int recv_count = send_counts[static_cast<size_t>(rank)];
  std::vector<int64_t> local_data(static_cast<size_t>(recv_count));

  if (world_size == 1) {
    local_data = GetInput();
  } else {
    MpiScattervInt64(GetInput(), send_counts, send_displs, local_data, recv_count, 0);
  }

  if (recv_count > 1) {
    ParallelRadixSort(local_data, 0, static_cast<size_t>(recv_count));
  }

  HierarchicalMerge(local_data, rank, world_size);

  BroadcastResult(local_data, rank, total);

  GetOutput() = std::move(local_data);

  return true;
}

}  // namespace leonova_a_radix_merge_sort
