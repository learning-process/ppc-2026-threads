#include "khruev_a_radix_sorting_int_bather_merge/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

#include "khruev_a_radix_sorting_int_bather_merge/common/include/common.hpp"

namespace khruev_a_radix_sorting_int_bather_merge {

void KhruevARadixSortingIntBatherMergeALL::CompareExchange(std::vector<int> &a, size_t i, size_t j) {
  if (a[i] > a[j]) {
    std::swap(a[i], a[j]);
  }
}

void KhruevARadixSortingIntBatherMergeALL::RadixSort(std::vector<int> &arr) {
  const int bits = 8;
  const int buckets = 1 << bits;
  const int mask = buckets - 1;
  const int passes = 32 / bits;

  std::vector<int> temp(arr.size());
  std::vector<int> *src = &arr;
  std::vector<int> *dst = &temp;

  for (int pass = 0; pass < passes; pass++) {
    std::vector<int> count(buckets, 0);
    int shift = pass * bits;

    for (int x : *src) {
      uint32_t ux = static_cast<uint32_t>(x) ^ 0x80000000U;
      uint32_t digit = (ux >> shift) & mask;
      count[digit]++;
    }

    for (int i = 1; i < buckets; i++) {
      count[i] += count[i - 1];
    }

    for (size_t i = src->size(); i-- > 0;) {
      uint32_t ux = static_cast<uint32_t>((*src)[i]) ^ 0x80000000U;
      uint32_t digit = (ux >> shift) & mask;
      (*dst)[--count[digit]] = (*src)[i];
    }

    std::swap(src, dst);
  }
}

void KhruevARadixSortingIntBatherMergeALL::OddEvenMerge(std::vector<int> &a, size_t n) {
  if (n < 2) {
    return;
  }

  size_t po = n / 2;

#pragma omp parallel for default(none) shared(a, po)
  for (size_t i = 0; i < po; ++i) {
    CompareExchange(a, i, i + po);
  }

  po >>= 1;

  for (; po > 0; po >>= 1) {
    size_t num_blocks = (n - (2 * po)) / (2 * po);

#pragma omp parallel for default(none) shared(a, po, num_blocks)
    for (size_t block_idx = 0; block_idx < num_blocks; ++block_idx) {
      size_t i = po + (block_idx * 2 * po);
      for (size_t j = 0; j < po; ++j) {
        CompareExchange(a, i + j, i + j + po);
      }
    }
  }
}

void KhruevARadixSortingIntBatherMergeALL::ScatterDataHelper(int rank, int active_procs, int chunk_size, size_t pow2,
                                                             std::vector<int> &local_data) {
  if (rank == 0) {
    std::vector<int> padded_data = GetInput();
    padded_data.resize(pow2, std::numeric_limits<int>::max());

    auto chunk_diff = static_cast<std::ptrdiff_t>(chunk_size);
    std::copy(padded_data.begin(), padded_data.begin() + chunk_diff, local_data.begin());

    for (int i = 1; i < active_procs; ++i) {
      auto offset = static_cast<std::ptrdiff_t>(i) * chunk_size;
      MPI_Send(padded_data.data() + offset, chunk_size, MPI_INT, i, 0, MPI_COMM_WORLD);
    }
  } else {
    MPI_Recv(local_data.data(), chunk_size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
}

void KhruevARadixSortingIntBatherMergeALL::LocalSortHelper(std::vector<int> &local_data) {
  size_t half_local = local_data.size() / 2;
  if (half_local > 0 && local_data.size() >= 2) {
    auto half_diff = static_cast<std::ptrdiff_t>(half_local);
    std::vector<int> left(local_data.begin(), local_data.begin() + half_diff);
    std::vector<int> right(local_data.begin() + half_diff, local_data.end());

#pragma omp parallel sections default(none) shared(left, right)
    {
#pragma omp section
      {
        RadixSort(left);
      }
#pragma omp section
      {
        RadixSort(right);
      }
    }

    std::ranges::copy(left, local_data.begin());
    std::ranges::copy(right, local_data.begin() + half_diff);
    OddEvenMerge(local_data, local_data.size());
  } else {
    RadixSort(local_data);
  }
}

void KhruevARadixSortingIntBatherMergeALL::TreeMergeHelper(int rank, int active_procs, int chunk_size,
                                                           std::vector<int> &local_data) {
  for (int step = 1; step < active_procs; step *= 2) {
    if (rank % (2 * step) == 0) {
      int sender = rank + step;
      int recv_count = chunk_size * step;
      std::vector<int> recv_data(recv_count);

      MPI_Recv(recv_data.data(), recv_count, MPI_INT, sender, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      local_data.insert(local_data.end(), recv_data.begin(), recv_data.end());
      OddEvenMerge(local_data, local_data.size());
    } else {
      int receiver = rank - step;
      int send_count = static_cast<int>(local_data.size());
      MPI_Send(local_data.data(), send_count, MPI_INT, receiver, 0, MPI_COMM_WORLD);
      break;
    }
  }
}

KhruevARadixSortingIntBatherMergeALL::KhruevARadixSortingIntBatherMergeALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool KhruevARadixSortingIntBatherMergeALL::ValidationImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int is_valid = 0;
  if (rank == 0) {
    is_valid = !GetInput().empty() ? 1 : 0;
  }
  MPI_Bcast(&is_valid, 1, MPI_INT, 0, MPI_COMM_WORLD);
  return is_valid == 1;
}

bool KhruevARadixSortingIntBatherMergeALL::PreProcessingImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    GetOutput().resize(GetInput().size());
  }
  return true;
}

bool KhruevARadixSortingIntBatherMergeALL::RunImpl() {
  int rank = 0;
  int num_procs = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  int original_size_int = 0;
  if (rank == 0) {
    original_size_int = static_cast<int>(GetInput().size());
  }

  MPI_Bcast(&original_size_int, 1, MPI_INT, 0, MPI_COMM_WORLD);
  auto original_size = static_cast<size_t>(original_size_int);

  if (original_size <= 1) {
    if (rank == 0) {
      GetOutput() = GetInput();
    } else {
      GetOutput().resize(original_size);
    }
    if (original_size == 1) {
      MPI_Bcast(GetOutput().data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    return true;
  }

  int active_procs = 1;
  while (active_procs * 2 <= num_procs) {
    active_procs *= 2;
  }

  std::vector<int> local_data;

  if (rank < active_procs) {
    size_t pow2 = 1;
    while (pow2 < original_size) {
      pow2 <<= 1;
    }
    while (std::cmp_less(pow2, active_procs)) {
      pow2 <<= 1;
    }

    int chunk_size = static_cast<int>(pow2 / active_procs);
    local_data.resize(chunk_size);

    ScatterDataHelper(rank, active_procs, chunk_size, pow2, local_data);
    LocalSortHelper(local_data);
    TreeMergeHelper(rank, active_procs, chunk_size, local_data);
  }

  local_data.resize(original_size);
  MPI_Bcast(local_data.data(), static_cast<int>(original_size), MPI_INT, 0, MPI_COMM_WORLD);
  GetOutput() = local_data;

  return true;
}

bool KhruevARadixSortingIntBatherMergeALL::PostProcessingImpl() {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    return GetOutput().size() == GetInput().size();
  }
  return true;
}

}  // namespace khruev_a_radix_sorting_int_bather_merge
