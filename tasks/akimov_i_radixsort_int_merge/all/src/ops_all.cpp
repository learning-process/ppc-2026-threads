#include "akimov_i_radixsort_int_merge/all/include/ops_all.hpp"

#include <mpi.h>
#include <oneapi/tbb/parallel_for.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <vector>

#include "akimov_i_radixsort_int_merge/common/include/common.hpp"
#include "util/include/util.hpp"

namespace akimov_i_radixsort_int_merge {

namespace {

void CountingSortStep(std::vector<int>::iterator in_begin, std::vector<int>::iterator in_end,
                      std::vector<int>::iterator out_begin, size_t byte_index) {
  size_t shift = byte_index * 8;
  std::vector<size_t> count(256, 0);

  for (auto it = in_begin; it != in_end; ++it) {
    auto raw_val = static_cast<uint32_t>(*it);
    uint32_t byte_val = (raw_val >> shift) & 0xFF;
    ++count[byte_val];
  }

  std::vector<size_t> prefix(256, 0);
  prefix[0] = 0;
  for (int i = 1; i < 256; ++i) {
    prefix[i] = prefix[i - 1] + count[i - 1];
  }

  for (auto it = in_begin; it != in_end; ++it) {
    auto raw_val = static_cast<uint32_t>(*it);
    uint32_t byte_val = (raw_val >> shift) & 0xFF;
    *(out_begin + static_cast<std::ptrdiff_t>(prefix[byte_val])) = *it;
    ++prefix[byte_val];
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
    tbb::parallel_for(0, num_blocks, [&, step](int i) {
      if (i % (2 * step) == 0 && i + step < num_blocks) {
        auto begin = arr.begin() + offsets[i];
        auto middle = arr.begin() + offsets[i + step];
        int end_idx = std::min(i + (2 * step), num_blocks);
        auto end = arr.begin() + offsets[end_idx];
        std::inplace_merge(begin, middle, end);
      }
    });
  }
}

void DistributeData(int rank, int world_size, int n, const std::vector<int> &input, std::vector<int> &local_data,
                    std::vector<int> &send_counts, std::vector<int> &send_displs) {
  send_counts.resize(world_size);
  send_displs.resize(world_size);
  int base = n / world_size;
  int rem = n % world_size;
  int offset = 0;
  for (int i = 0; i < world_size; ++i) {
    send_counts[i] = base + (i < rem ? 1 : 0);
    send_displs[i] = offset;
    offset += send_counts[i];
  }

  int local_size = send_counts[rank];
  local_data.resize(local_size);

  if (ppc::util::IsUnderMpirun()) {
    MPI_Scatterv(rank == 0 ? input.data() : nullptr, send_counts.data(), send_displs.data(), MPI_INT, local_data.data(),
                 local_size, MPI_INT, 0, MPI_COMM_WORLD);
  } else {
    local_data = input;
  }
}

void GatherData(int rank, int n, const std::vector<int> &local_data, const std::vector<int> &send_counts,
                const std::vector<int> &send_displs, std::vector<int> &output) {
  if (ppc::util::IsUnderMpirun()) {
    if (rank == 0) {
      output.resize(n);
    }
    MPI_Gatherv(local_data.data(), static_cast<int>(local_data.size()), MPI_INT, rank == 0 ? output.data() : nullptr,
                send_counts.data(), send_displs.data(), MPI_INT, 0, MPI_COMM_WORLD);
  } else {
    output = local_data;
  }
}

void InvertSignBit(std::vector<int> &data, int32_t mask) {
  for (int &val : data) {
    val ^= mask;
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
  int rank = 0;
  if (ppc::util::IsUnderMpirun()) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  }
  if (rank == 0) {
    GetOutput() = GetInput();
  }
  return true;
}

bool AkimovIRadixSortIntMergeALL::RunImpl() {
  const bool is_mpi = ppc::util::IsUnderMpirun();
  int rank = 0;
  int world_size = 1;

  if (is_mpi) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  }

  int n = 0;
  if (rank == 0) {
    n = static_cast<int>(GetOutput().size());
  }

  if (is_mpi) {
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }

  if (n == 0) {
    return true;
  }

  std::vector<int> send_counts;
  std::vector<int> send_displs;
  std::vector<int> local_data;
  DistributeData(rank, world_size, n, GetOutput(), local_data, send_counts, send_displs);

  constexpr int32_t kSignMask = INT32_MIN;
  InvertSignBit(local_data, kSignMask);
  RadixSortLocal(local_data.begin(), local_data.end());
  InvertSignBit(local_data, kSignMask);

  GatherData(rank, n, local_data, send_counts, send_displs, GetOutput());

  if (rank == 0 && world_size > 1) {
    std::vector<int> offsets(world_size + 1, 0);
    for (int i = 0; i < world_size; ++i) {
      offsets[i + 1] = offsets[i] + send_counts[i];
    }
    ParallelMerge(GetOutput(), offsets, world_size);
  }

  if (is_mpi) {
    int output_size = static_cast<int>(GetOutput().size());
    MPI_Bcast(&output_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0) {
      GetOutput().resize(output_size);
    }
    MPI_Bcast(GetOutput().data(), output_size, MPI_INT, 0, MPI_COMM_WORLD);
  }

  return true;
}

bool AkimovIRadixSortIntMergeALL::PostProcessingImpl() {
  return true;
}

}  // namespace akimov_i_radixsort_int_merge
