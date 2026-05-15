#include "frolova_s_radix_sort_double/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "frolova_s_radix_sort_double/common/include/common.hpp"

namespace frolova_s_radix_sort_double {

namespace {

void ParallelRadixSortOpenMP(std::vector<double> &data) {
  const std::size_t n = data.size();
  if (n < 2) {
    return;
  }

  constexpr int kRadix = 256;
  constexpr int kNumBits = 8;
  constexpr int kNumPasses = sizeof(std::uint64_t);

  std::vector<double> temp(n);

  for (int pass = 0; pass < kNumPasses; ++pass) {
    std::vector<std::vector<int>> local_counts(omp_get_max_threads(), std::vector<int>(kRadix, 0));

#pragma omp parallel
    {
      int tid = omp_get_thread_num();
      std::vector<int> &count = local_counts[tid];

#pragma omp for schedule(static)
      for (std::size_t i = 0; i < n; ++i) {
        auto bits = std::bit_cast<std::uint64_t>(data[i]);
        int byte = static_cast<int>((bits >> (pass * kNumBits)) & 0xFF);
        ++count[byte];
      }
    }

    std::vector<int> total_count(kRadix, 0);
    for (const auto &local : local_counts) {
      for (int i = 0; i < kRadix; ++i) {
        total_count[i] += local[i];
      }
    }

    std::vector<int> offsets(kRadix, 0);
    int sum = 0;
    for (int i = 0; i < kRadix; ++i) {
      offsets[i] = sum;
      sum += total_count[i];
    }

    std::vector<int> current_offsets = offsets;
#pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < n; ++i) {
      auto bits = std::bit_cast<std::uint64_t>(data[i]);
      int byte = static_cast<int>((bits >> (pass * kNumBits)) & 0xFF);
      int pos = 0;
#pragma omp atomic capture
      pos = current_offsets[byte]++;
      temp[pos] = data[i];
    }

    data.swap(temp);
  }
}

void FixNegativeOrderOpenMP(std::vector<double> &data) {
  std::size_t n = data.size();
  std::size_t first_non_negative = 0;
  while (first_non_negative < n && (std::bit_cast<std::uint64_t>(data[first_non_negative]) >> 63) != 0) {
    ++first_non_negative;
  }
  std::size_t neg_count = first_non_negative;
#pragma omp parallel for
  for (std::size_t i = 0; i < neg_count / 2; ++i) {
    std::swap(data[i], data[neg_count - 1 - i]);
  }
}

}  // namespace

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
  int mpi_initialized_flag = 0;
  MPI_Initialized(&mpi_initialized_flag);

  if (!mpi_initialized_flag) {
    // Запуск без MPI — только OpenMP
    std::vector<double> working = GetInput();
    ParallelRadixSortOpenMP(working);
    FixNegativeOrderOpenMP(working);
    GetOutput() = std::move(working);
    return true;
  }

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const InType &full_input = GetInput();
  std::size_t total_n = full_input.size();

  std::vector<int> sendcounts(size, 0);
  std::vector<int> displs(size, 0);

  std::size_t base = total_n / size;
  int remainder = total_n % size;

  for (int i = 0; i < size; ++i) {
    sendcounts[i] = static_cast<int>(base + (i < remainder ? 1 : 0));
    displs[i] = (i == 0) ? 0 : displs[i - 1] + sendcounts[i - 1];
  }

  std::vector<double> local_data(sendcounts[rank]);
  MPI_Scatterv(full_input.data(), sendcounts.data(), displs.data(), MPI_DOUBLE, local_data.data(), sendcounts[rank],
               MPI_DOUBLE, 0, MPI_COMM_WORLD);

  ParallelRadixSortOpenMP(local_data);
  FixNegativeOrderOpenMP(local_data);

  if (rank == 0) {
    std::vector<double> gathered(total_n);
    MPI_Gatherv(local_data.data(), local_data.size(), MPI_DOUBLE, gathered.data(), sendcounts.data(), displs.data(),
                MPI_DOUBLE, 0, MPI_COMM_WORLD);

    std::vector<std::size_t> indices(size, 0);
    OutType result;
    result.reserve(total_n);

    while (result.size() < total_n) {
      int best_proc = -1;
      double min_val = 0.0;
      for (int p = 0; p < size; ++p) {
        if (indices[p] < static_cast<std::size_t>(sendcounts[p])) {
          double val = gathered[displs[p] + indices[p]];
          if (best_proc == -1 || val < min_val) {
            min_val = val;
            best_proc = p;
          }
        }
      }
      if (best_proc == -1) {
        break;
      }
      result.push_back(min_val);
      ++indices[best_proc];
    }

    GetOutput() = std::move(result);
  } else {
    MPI_Gatherv(local_data.data(), local_data.size(), MPI_DOUBLE, nullptr, nullptr, nullptr, MPI_DOUBLE, 0,
                MPI_COMM_WORLD);
    GetOutput().resize(total_n);
  }

  MPI_Bcast(GetOutput().data(), total_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  return true;
}

bool FrolovaSRadixSortDoubleALL::PostProcessingImpl() {
  return true;
}

}  // namespace frolova_s_radix_sort_double
