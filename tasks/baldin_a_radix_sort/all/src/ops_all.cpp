#include "baldin_a_radix_sort/all/include/ops_all.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <vector>

#include "baldin_a_radix_sort/common/include/common.hpp"
#include "util/include/util.hpp"

namespace baldin_a_radix_sort {  // comment for ci rerun

BaldinARadixSortALL::BaldinARadixSortALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool BaldinARadixSortALL::ValidationImpl() {
  return true;
}

bool BaldinARadixSortALL::PreProcessingImpl() {
  int rank = 0;

  if (ppc::util::IsUnderMpirun()) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  }

  if (rank == 0) {
    GetOutput() = GetInput();
  }

  return true;
}

namespace {

void CountingSortStep(std::vector<int>::iterator in_begin, std::vector<int>::iterator in_end,
                      std::vector<int>::iterator out_begin, size_t byte_index) {
  size_t shift = byte_index * 8;
  std::array<size_t, 256> count = {0};

  for (auto it = in_begin; it != in_end; it++) {
    auto raw_val = static_cast<unsigned int>(*it);
    unsigned int byte_val = (raw_val >> shift) & 0xFF;

    if (byte_index == sizeof(int) - 1) {
      byte_val ^= 128;
    }
    count.at(byte_val)++;
  }

  std::array<size_t, 256> prefix{};
  prefix[0] = 0;
  for (int i = 1; i < 256; i++) {
    prefix.at(i) = prefix.at(i - 1) + count.at(i - 1);
  }

  for (auto it = in_begin; it != in_end; it++) {
    auto raw_val = static_cast<unsigned int>(*it);
    unsigned int byte_val = (raw_val >> shift) & 0xFF;

    if (byte_index == sizeof(int) - 1) {
      byte_val ^= 128;
    }

    *(out_begin + static_cast<int64_t>(prefix.at(byte_val))) = *it;
    prefix.at(byte_val)++;
  }
}

void RadixSortLocal(std::vector<int>::iterator begin, std::vector<int>::iterator end) {
  size_t n = std::distance(begin, end);
  if (n < 2) {
    return;
  }

  std::vector<int> temp(n);

  for (size_t i = 0; i < sizeof(int); i++) {
    size_t shift = i;

    if (i % 2 == 0) {
      CountingSortStep(begin, end, temp.begin(), shift);
    } else {
      CountingSortStep(temp.begin(), temp.end(), begin, shift);
    }
  }
}

}  // namespace

void BaldinARadixSortALL::DistributeData(int rank, int size, int &n) {
  bool is_mpi = ppc::util::IsUnderMpirun();

  if (is_mpi) {
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }

  if (n == 0) {
    return;
  }

  counts_.resize(size);
  displs_.resize(size);
  int chunk = n / size;
  int rem = n % size;
  int curr = 0;
  for (int i = 0; i < size; i++) {
    counts_[i] = chunk + (i < rem ? 1 : 0);
    displs_[i] = curr;
    curr += counts_[i];
  }

  local_data_.resize(counts_[rank]);

  if (is_mpi) {
    MPI_Scatterv(rank == 0 ? GetOutput().data() : nullptr, counts_.data(), displs_.data(), MPI_INT, local_data_.data(),
                 counts_[rank], MPI_INT, 0, MPI_COMM_WORLD);
  } else {
    local_data_ = GetOutput();
  }
}

void BaldinARadixSortALL::GatherData(int rank) {
  bool is_mpi = ppc::util::IsUnderMpirun();

  if (is_mpi) {
    MPI_Gatherv(local_data_.data(), counts_[rank], MPI_INT, rank == 0 ? GetOutput().data() : nullptr, counts_.data(),
                displs_.data(), MPI_INT, 0, MPI_COMM_WORLD);
  } else {
    GetOutput() = local_data_;
  }
}

void BaldinARadixSortALL::LocalSort(int num_threads) {
  int local_n = static_cast<int>(local_data_.size());
  if (local_n == 0) {
    return;
  }

  std::vector<int> local_offsets(num_threads + 1);
  int l_chunk = local_n / num_threads;
  int l_rem = local_n % num_threads;
  int l_curr = 0;
  for (int i = 0; i < num_threads; i++) {
    local_offsets[i] = l_curr;
    l_curr += l_chunk + (i < l_rem ? 1 : 0);
  }
  local_offsets[num_threads] = local_n;

  auto &local_data_ref = local_data_;

#pragma omp parallel num_threads(num_threads) default(none) shared(num_threads, local_offsets, local_data_ref)
  {
    int tid = omp_get_thread_num();
    auto begin = local_data_ref.begin() + local_offsets[tid];
    auto end = local_data_ref.begin() + local_offsets[tid + 1];
    RadixSortLocal(begin, end);
  }

  for (int step = 1; step < num_threads; step *= 2) {
#pragma omp parallel for num_threads(num_threads) default(none) shared(step, num_threads, local_offsets, local_data_ref)
    for (int i = 0; i < num_threads; i += 2 * step) {
      if (i + step < num_threads) {
        auto begin = local_data_ref.begin() + local_offsets[i];
        auto middle = local_data_ref.begin() + local_offsets[i + step];
        int end_idx = std::min(i + (2 * step), num_threads);
        auto end = local_data_ref.begin() + local_offsets[end_idx];

        std::inplace_merge(begin, middle, end);
      }
    }
  }
}

void BaldinARadixSortALL::GlobalMerge(int num_threads, int size, int n) {
  auto &out = GetOutput();

  std::vector<int> global_offsets(size + 1);
  for (int i = 0; i < size; i++) {
    global_offsets[i] = displs_[i];
  }
  global_offsets[size] = n;

  for (int step = 1; step < size; step *= 2) {
#pragma omp parallel for num_threads(num_threads) default(none) shared(step, size, global_offsets, out)
    for (int i = 0; i < size; i += 2 * step) {
      if (i + step < size) {
        auto begin = out.begin() + global_offsets[i];
        auto middle = out.begin() + global_offsets[i + step];
        int end_idx = std::min(i + (2 * step), size);
        auto end = out.begin() + global_offsets[end_idx];

        std::inplace_merge(begin, middle, end);
      }
    }
  }
}

bool BaldinARadixSortALL::RunImpl() {
  int rank = 0;
  int size = 1;
  bool is_mpi = ppc::util::IsUnderMpirun();

  if (is_mpi) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
  }

  int n = 0;
  if (rank == 0) {
    n = static_cast<int>(GetOutput().size());
  }

  DistributeData(rank, size, n);

  if (n == 0) {
    return true;
  }

  int num_threads = ppc::util::GetNumThreads();

  LocalSort(num_threads);

  GatherData(rank);

  if (rank == 0 && size > 1) {
    GlobalMerge(num_threads, size, n);
  }

  return true;
}

bool BaldinARadixSortALL::PostProcessingImpl() {
  return true;
}

}  // namespace baldin_a_radix_sort
