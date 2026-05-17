#include "akimov_i_radixsort_int_merge/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <utility>
#include <vector>

#include "akimov_i_radixsort_int_merge/common/include/common.hpp"

namespace akimov_i_radixsort_int_merge {

namespace {

void CountingSortStep(std::vector<int>::iterator in_begin, std::vector<int>::iterator in_end,
                      std::vector<int>::iterator out_begin, size_t byte_index) {
  size_t shift = byte_index * 8;
  std::array<size_t, 256> count = {0};

  for (auto it = in_begin; it != in_end; ++it) {
    auto raw_val = static_cast<unsigned int>(*it);
    unsigned int byte_val = (raw_val >> shift) & 0xFF;
    count.at(byte_val)++;
  }

  std::array<size_t, 256> prefix{};
  prefix[0] = 0;
  for (int i = 1; i < 256; ++i) {
    prefix.at(i) = prefix.at(i - 1) + count.at(i - 1);
  }

  for (auto it = in_begin; it != in_end; ++it) {
    auto raw_val = static_cast<unsigned int>(*it);
    unsigned int byte_val = (raw_val >> shift) & 0xFF;
    *(out_begin + static_cast<std::ptrdiff_t>(prefix.at(byte_val))) = *it;
    prefix.at(byte_val)++;
  }
}

void RadixSortLocal(std::vector<int>::iterator begin, std::vector<int>::iterator end) {
  size_t n = std::distance(begin, end);
  if (n < 2) {
    return;
  }

  std::vector<int> temp(n);
  for (size_t i = 0; i < sizeof(int); ++i) {
    if (i % 2 == 0) {
      CountingSortStep(begin, end, temp.begin(), i);
    } else {
      CountingSortStep(temp.begin(), temp.end(), begin, i);
    }
  }
}

void ParallelMerge(std::vector<int> &arr, const std::vector<int> &offsets, int num_blocks) {
  for (int step = 1; step < num_blocks; step *= 2) {
    for (int i = 0; i < num_blocks; i += 2 * step) {
      if (i + step < num_blocks) {
        auto begin = arr.begin() + offsets[i];
        auto middle = arr.begin() + offsets[i + step];
        int end_idx = std::min(i + (2 * step), num_blocks);
        auto end = arr.begin() + offsets[end_idx];
        std::inplace_merge(begin, middle, end);
      }
    }
  }
}

}  // namespace

AkimovIRadixSortIntMergeALL::AkimovIRadixSortIntMergeALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType();
}

bool AkimovIRadixSortIntMergeALL::ValidationImpl() {
  return !GetInput().empty();
}

bool AkimovIRadixSortIntMergeALL::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool AkimovIRadixSortIntMergeALL::RunImpl() {
  auto &arr = GetOutput();
  int n = static_cast<int>(arr.size());
  if (n == 0) {
    return true;
  }

  int rank = -1;
  int world_size = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  std::vector<int> send_counts(world_size);
  std::vector<int> send_displs(world_size);
  int base = n / world_size;
  int rem = n % world_size;
  int offset = 0;
  for (int i = 0; i < world_size; ++i) {
    send_counts[i] = base + (i < rem ? 1 : 0);
    send_displs[i] = offset;
    offset += send_counts[i];
  }

  int local_size = send_counts[rank];
  std::vector<int> local_data(local_size);

  MPI_Scatterv(arr.data(), send_counts.data(), send_displs.data(), MPI_INT, local_data.data(), local_size, MPI_INT, 0,
               MPI_COMM_WORLD);

  constexpr int32_t kSignMask = INT32_MIN;
  for (int i = 0; i < local_size; ++i) {
    local_data[i] ^= kSignMask;
  }

  RadixSortLocal(local_data.begin(), local_data.end());

  for (int i = 0; i < local_size; ++i) {
    local_data[i] ^= kSignMask;
  }

  std::vector<int> global_data;
  if (rank == 0) {
    global_data.resize(n);
  }
  MPI_Gatherv(local_data.data(), local_size, MPI_INT, global_data.data(), send_counts.data(), send_displs.data(),
              MPI_INT, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    std::vector<int> offsets(world_size + 1, 0);
    for (int i = 0; i < world_size; ++i) {
      offsets[i + 1] = offsets[i] + send_counts[i];
    }
    ParallelMerge(global_data, offsets, world_size);
    GetOutput() = std::move(global_data);
  } else {
    GetOutput().clear();
  }

  int output_size = static_cast<int>(GetOutput().size());
  MPI_Bcast(&output_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank != 0) {
    GetOutput().resize(output_size);
  }
  MPI_Bcast(GetOutput().data(), output_size, MPI_INT, 0, MPI_COMM_WORLD);

  return true;
}

bool AkimovIRadixSortIntMergeALL::PostProcessingImpl() {
  return true;
}

}  // namespace akimov_i_radixsort_int_merge
