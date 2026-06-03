#include "dorofeev_i_bitwise_sort_double_eo_batcher_merge/all/include/ops_all.hpp"

#include <mpi.h>
#include <tbb/tbb.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

#include "dorofeev_i_bitwise_sort_double_eo_batcher_merge/common/include/common.hpp"
#include "util/include/util.hpp"

namespace dorofeev_i_bitwise_sort_double_eo_batcher_merge {

namespace {

uint64_t DoubleToUint(double d) {
  uint64_t u = 0;
  std::memcpy(&u, &d, sizeof(double));
  if ((u & 0x8000000000000000ULL) != 0) {
    u = ~u;
  } else {
    u |= 0x8000000000000000ULL;
  }
  return u;
}

double UintToDouble(uint64_t u) {
  if ((u & 0x8000000000000000ULL) != 0) {
    u &= ~0x8000000000000000ULL;
  } else {
    u = ~u;
  }
  double d = 0.0;
  std::memcpy(&d, &u, sizeof(double));
  return d;
}

void RadixSortDouble(std::vector<double> &arr) {
  if (arr.empty()) {
    return;
  }

  std::vector<uint64_t> uarr(arr.size());
  for (size_t i = 0; i < arr.size(); ++i) {
    uarr[i] = DoubleToUint(arr[i]);
  }

  std::vector<uint64_t> temp(uarr.size());
  for (size_t byte = 0; byte < 8; ++byte) {
    std::vector<int> count(256, 0);
    for (uint64_t val : uarr) {
      count[(val >> (byte * 8)) & 0xFF]++;
    }
    for (size_t i = 1; i < 256; ++i) {
      count[i] += count[i - 1];
    }
    for (int i = static_cast<int>(uarr.size()) - 1; i >= 0; --i) {
      temp[--count[(uarr[i] >> (byte * 8)) & 0xFF]] = uarr[i];
    }
    uarr = temp;
  }

  for (size_t i = 0; i < arr.size(); ++i) {
    arr[i] = UintToDouble(uarr[i]);
  }
}

void CompareExchangeBlocks(double *arr, size_t i, size_t step) {
  for (size_t k = 0; k < step; ++k) {
    if (arr[i + k] > arr[i + k + step]) {
      std::swap(arr[i + k], arr[i + k + step]);
    }
  }
}

void OddEvenMergeIterative(double *arr, size_t start, size_t n) {
  if (n <= 1) {
    return;
  }
  size_t step = n / 2;
  CompareExchangeBlocks(arr, start, step);
  step /= 2;
  for (; step > 0; step /= 2) {
    for (size_t i = step; i < n - step; i += step * 2) {
      CompareExchangeBlocks(arr, start + i, step);
    }
  }
}

void ProcessChunkTBB(double *raw_data, int chunk_idx, size_t chunk_size) {
  size_t start_idx = static_cast<size_t>(chunk_idx) * chunk_size;
  std::vector<double> local_arr(chunk_size);
  for (size_t j = 0; j < chunk_size; ++j) {
    local_arr[j] = raw_data[start_idx + j];
  }
  RadixSortDouble(local_arr);
  for (size_t j = 0; j < chunk_size; ++j) {
    raw_data[start_idx + j] = local_arr[j];
  }
}

void ExecuteTBBSortLocal(double *raw_data, size_t total_size, size_t chunk_size, int num_chunks_int) {
  int num_threads = ppc::util::GetNumThreads();
  tbb::task_arena arena(num_threads > 0 ? num_threads : 1);
  arena.execute([&] {
    tbb::parallel_for(tbb::blocked_range<int>(0, num_chunks_int),
                      [raw_data, chunk_size](const tbb::blocked_range<int> &r) {
      for (int i = r.begin(); i != r.end(); ++i) {
        ProcessChunkTBB(raw_data, i, chunk_size);
      }
    });
    for (size_t size = chunk_size; size < total_size; size *= 2) {
      int merges_count = static_cast<int>(total_size / (size * 2));
      tbb::parallel_for(tbb::blocked_range<int>(0, merges_count), [raw_data, size](const tbb::blocked_range<int> &r) {
        for (int i = r.begin(); i != r.end(); ++i) {
          OddEvenMergeIterative(raw_data, static_cast<size_t>(i) * 2 * size, 2 * size);
        }
      });
    }
  });
}

void PerformLocalTBBSort(double *data, size_t total_size) {
  if (total_size == 0) {
    return;
  }
  int num_threads = ppc::util::GetNumThreads();
  if (num_threads <= 0) {
    num_threads = 1;
  }
  size_t num_chunks = 1;
  while (num_chunks * 2 <= static_cast<size_t>(num_threads) && num_chunks * 2 <= total_size) {
    num_chunks *= 2;
  }
  ExecuteTBBSortLocal(data, total_size, total_size / num_chunks, static_cast<int>(num_chunks));
}

void PerformFinalMergeTBB(double *raw_data, size_t local_chunk_size, size_t pow2) {
  int num_threads = ppc::util::GetNumThreads();
  tbb::task_arena arena(num_threads > 0 ? num_threads : 1);
  arena.execute([&] {
    for (size_t size = local_chunk_size; size < pow2; size *= 2) {
      int merges_count = static_cast<int>(pow2 / (size * 2));
      tbb::parallel_for(tbb::blocked_range<int>(0, merges_count), [raw_data, size](const tbb::blocked_range<int> &r) {
        for (int i = r.begin(); i != r.end(); ++i) {
          OddEvenMergeIterative(raw_data, static_cast<size_t>(i) * 2 * size, 2 * size);
        }
      });
    }
  });
}

void ExecuteMPIHybridSort(std::vector<double> &local_data, size_t pow2, int world_rank, int active_procs) {
  size_t local_chunk_size = pow2 / active_procs;
  int chunk_size_int = static_cast<int>(local_chunk_size);
  size_t buffer_size = (world_rank < active_procs) ? local_chunk_size : 0;
  std::vector<double> mpi_buffer(buffer_size, 0.0);

  if (world_rank == 0) {
    MPI_Scatter(local_data.data(), chunk_size_int, MPI_DOUBLE, mpi_buffer.data(), chunk_size_int, MPI_DOUBLE, 0,
                MPI_COMM_WORLD);
  } else if (world_rank < active_procs) {
    MPI_Scatter(nullptr, 0, MPI_DATATYPE_NULL, mpi_buffer.data(), chunk_size_int, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }

  if (world_rank < active_procs) {
    PerformLocalTBBSort(mpi_buffer.data(), local_chunk_size);
  }

  if (world_rank == 0) {
    MPI_Gather(mpi_buffer.data(), chunk_size_int, MPI_DOUBLE, local_data.data(), chunk_size_int, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
    PerformFinalMergeTBB(local_data.data(), local_chunk_size, pow2);
  } else if (world_rank < active_procs) {
    MPI_Gather(mpi_buffer.data(), chunk_size_int, MPI_DOUBLE, nullptr, 0, MPI_DATATYPE_NULL, 0, MPI_COMM_WORLD);
  }
}

void HandleGTestLocalSort(std::vector<double> &local_data, int world_rank) {
  if (world_rank == 0 || local_data.empty()) {
    return;
  }
  size_t my_orig = local_data.size();
  size_t my_pow2 = 1;
  while (my_pow2 < my_orig) {
    my_pow2 *= 2;
  }
  if (my_pow2 > my_orig) {
    local_data.resize(my_pow2, std::numeric_limits<double>::max());
  }

  PerformLocalTBBSort(local_data.data(), my_pow2);

  if (my_pow2 > my_orig) {
    local_data.resize(my_orig);
  }
}

size_t CalculatePaddedSize(size_t original_size) {
  if (original_size == 0) {
    return 1;
  }
  size_t pow2 = 1;
  while (pow2 < original_size) {
    pow2 *= 2;
  }
  return pow2;
}

}  // namespace

DorofeevIBitwiseSortDoubleEOBatcherMergeALL::DorofeevIBitwiseSortDoubleEOBatcherMergeALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool DorofeevIBitwiseSortDoubleEOBatcherMergeALL::ValidationImpl() {
  return true;
}

bool DorofeevIBitwiseSortDoubleEOBatcherMergeALL::PreProcessingImpl() {
  local_data_ = GetInput();
  return true;
}

bool DorofeevIBitwiseSortDoubleEOBatcherMergeALL::RunImpl() {
  int world_size = 0;
  int world_rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  bool is_pow2_procs = (world_size > 0) && ((world_size & (world_size - 1)) == 0);
  int active_procs = is_pow2_procs ? world_size : 1;

  size_t rank0_original_size = 0;
  size_t pow2 = 1;

  if (world_rank == 0) {
    rank0_original_size = local_data_.size();
    if (rank0_original_size == 0) {
      active_procs = 1;
    }
    pow2 = CalculatePaddedSize(rank0_original_size);
    if (pow2 > rank0_original_size) {
      local_data_.resize(pow2, std::numeric_limits<double>::max());
    }
  }

  MPI_Bcast(&active_procs, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&pow2, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);

  if (active_procs == 1) {
    if (world_rank == 0 && pow2 > 0) {
      PerformLocalTBBSort(local_data_.data(), pow2);
    }
  } else {
    ExecuteMPIHybridSort(local_data_, pow2, world_rank, active_procs);
  }

  if (world_rank == 0 && pow2 > rank0_original_size) {
    local_data_.resize(rank0_original_size);
  }

  HandleGTestLocalSort(local_data_, world_rank);

  return true;
}

bool DorofeevIBitwiseSortDoubleEOBatcherMergeALL::PostProcessingImpl() {
  GetOutput() = local_data_;
  return true;
}

}  // namespace dorofeev_i_bitwise_sort_double_eo_batcher_merge
