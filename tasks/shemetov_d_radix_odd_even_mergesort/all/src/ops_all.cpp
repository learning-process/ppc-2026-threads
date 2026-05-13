#include "shemetov_d_radix_odd_even_mergesort/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <climits>
#include <cstddef>
#include <utility>
#include <vector>

#include "shemetov_d_radix_odd_even_mergesort/common/include/common.hpp"

namespace shemetov_d_radix_odd_even_mergesort {

ShemetovDRadixOddEvenMergeSortALL::ShemetovDRadixOddEvenMergeSortALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

void ShemetovDRadixOddEvenMergeSortALL::ScatterData(std::vector<int> &global_array, std::vector<int> &local_array,
                                                    size_t chunk, int rank, int ranks, int mpi_init) {
  if (mpi_init != 0 && ranks > 1) {
    if (rank == 0) {
      std::ranges::copy(global_array.begin(), global_array.begin() + static_cast<int>(chunk), local_array.begin());

      for (int rank_index = 1; rank_index < ranks; rank_index += 1) {
        MPI_Send(global_array.data() + (rank_index * chunk), static_cast<int>(chunk), MPI_INT, rank_index, 0,
                 MPI_COMM_WORLD);
      }
    } else {
      MPI_Recv(local_array.data(), static_cast<int>(chunk), MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
  } else if (rank == 0) {
    local_array = global_array;
  }
}

void ShemetovDRadixOddEvenMergeSortALL::MPISort(std::vector<int> &local_array, size_t chunk) {
  size_t threads = omp_get_max_threads();

  size_t limit = 1;
  while (limit * 2 <= std::min(threads, chunk)) {
    limit *= 2;
  }

  size_t chunk_size = chunk / limit;

#pragma omp parallel num_threads(limit) default(none) shared(local_array, chunk, chunk_size, limit)
  {
    size_t thread_num = omp_get_thread_num();
    size_t left = thread_num * chunk_size;
    size_t right = left + chunk_size - 1;

    std::vector<int> buffer;
    std::vector<int> position;

    RadixSort(local_array, left, right, buffer, position);
#pragma omp barrier
    for (size_t segment = chunk_size * 2; segment <= chunk; segment *= 2) {
#pragma omp for
      for (size_t i = 0; i < chunk; i += segment) {
        OddEvenMerge(local_array, i, segment);
      }
#pragma omp barrier
    }
  }
}

void ShemetovDRadixOddEvenMergeSortALL::MPIMerge(size_t chunk) {
  std::vector<int> &ref_array = array_;
  size_t power = power_;

  for (size_t segment = chunk * 2; segment <= power; segment *= 2) {
#pragma omp parallel for default(none) shared(ref_array, power, segment)
    for (size_t i = 0; i < power; i += segment) {
      OddEvenMerge(ref_array, i, segment);
    }
  }
}

void ShemetovDRadixOddEvenMergeSortALL::GatherData(std::vector<int> &global_array, std::vector<int> &local_array,
                                                   size_t chunk, int rank, int ranks, int mpi_init) {
  if (mpi_init != 0 && ranks > 1) {
    if (rank == 0) {
      std::ranges::copy(local_array, global_array.begin());

      for (int rank_index = 1; rank_index < ranks; rank_index += 1) {
        MPI_Recv(global_array.data() + (rank_index * chunk), static_cast<int>(chunk), MPI_INT, rank_index, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
    } else {
      MPI_Send(local_array.data(), static_cast<int>(chunk), MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
  } else if (rank == 0) {
    global_array = local_array;
  }
}

void ShemetovDRadixOddEvenMergeSortALL::RadixSort(std::vector<int> &array, size_t left, size_t right,
                                                  std::vector<int> &buffer, std::vector<int> &position) {
  if (left >= right) {
    return;
  }

  int maximum =
      *std::ranges::max_element(array.begin() + static_cast<int>(left), array.begin() + static_cast<int>(right) + 1);

  size_t segment = right - left + 1;

  buffer.resize(segment);
  for (size_t merge_shift = 0; merge_shift < 32; merge_shift += 8) {
    position.assign(256, 0);

    for (size_t i = left; i <= right; i += 1) {
      int apply_bitmask = (array[i] >> merge_shift) & 0xFF;

      position[apply_bitmask] += 1;
    }

    for (size_t i = 1; i < 256; i += 1) {
      position[i] += position[i - 1];
    }

    for (size_t i = segment; i > 0; i -= 1) {
      size_t current_index = left + i - 1;
      int apply_bitmask = (array[current_index] >> merge_shift) & 0xFF;

      buffer[position[apply_bitmask] -= 1] = array[current_index];
    }

    for (size_t i = 0; i < segment; i += 1) {
      array[left + i] = buffer[i];
    }

    if ((maximum >> merge_shift) < 256) {
      break;
    }
  }
}

void ShemetovDRadixOddEvenMergeSortALL::OddEvenMerge(std::vector<int> &array, size_t start_offset, size_t segment) {
  if (segment <= 1) {
    return;
  }

  size_t padding = segment / 2;
  for (size_t i = 0; i < padding; i += 1) {
    if (array[start_offset + i] > array[start_offset + padding + i]) {
      std::swap(array[start_offset + i], array[start_offset + padding + i]);
    }
  }

  for (padding = segment / 4; padding > 0; padding /= 2) {
    size_t step = padding * 2;

    for (size_t start_position = padding; start_position < segment - padding; start_position += step) {
      for (size_t i = 0; i < padding; i += 1) {
        if (array[start_offset + start_position + i] > array[start_offset + start_position + i + padding]) {
          std::swap(array[start_offset + start_position + i], array[start_offset + start_position + i + padding]);
        }
      }
    }
  }
}

bool ShemetovDRadixOddEvenMergeSortALL::ValidationImpl() {
  const auto &[size, array] = GetInput();
  return size > 0 && static_cast<size_t>(size) == array.size();
}

bool ShemetovDRadixOddEvenMergeSortALL::PreProcessingImpl() {
  const auto &[size, array] = GetInput();

  if (size <= 0) {
    return true;
  }

  array_ = array;
  offset_ = *std::ranges::min_element(array_.begin(), array_.end());
  size_ = static_cast<size_t>(size);
  power_ = 1;

  while (power_ < size_) {
    power_ *= 2;
  }

  for (size_t i = 0; i < size_; i += 1) {
    array_[i] -= offset_;
  }

  if (power_ > size_) {
    array_.resize(power_, INT_MAX);
  }

  return true;
}

bool ShemetovDRadixOddEvenMergeSortALL::RunImpl() {
  if (size_ == 0 || power_ <= 1) {
    return true;
  }

  int mpi_init = 0;
  MPI_Initialized(&mpi_init);

  int rank = 0;
  int num_procs = 1;

  if (mpi_init != 0) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  }

  int ranks = 1;
  while (ranks * 2 <= num_procs && std::cmp_less_equal(ranks * 2, power_)) {
    ranks *= 2;
  }

  if (rank < ranks) {
    size_t chunk = power_ / ranks;
    std::vector<int> segment_array(chunk);

    ScatterData(array_, segment_array, chunk, rank, ranks, mpi_init);

    MPISort(segment_array, chunk);

    GatherData(array_, segment_array, chunk, rank, ranks, mpi_init);

    if (rank == 0) {
      MPIMerge(chunk);
    }
  }

  if (mpi_init != 0) {
    MPI_Bcast(array_.data(), static_cast<int>(power_), MPI_INT, 0, MPI_COMM_WORLD);
  }

  return true;
}

bool ShemetovDRadixOddEvenMergeSortALL::PostProcessingImpl() {
  if (size_ == 0) {
    return true;
  }

  array_.resize(size_);

  for (size_t i = 0; i < size_; i += 1) {
    array_[i] += offset_;
  }

  if (!std::ranges::is_sorted(array_.begin(), array_.end())) {
    return false;
  }

  GetOutput() = array_;
  return true;
}

}  // namespace shemetov_d_radix_odd_even_mergesort
