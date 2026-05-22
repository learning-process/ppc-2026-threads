#include "gonozov_l_bitwise_sorting_double_Batcher_merge/all/include/ops_all.hpp"

#include <mpi.h>
#include <tbb/tbb.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

#include "gonozov_l_bitwise_sorting_double_Batcher_merge/common/include/common.hpp"
#include "util/include/util.hpp"

namespace gonozov_l_bitwise_sorting_double_batcher_merge {

namespace {

uint64_t DoubleToSortableInt(double d) {
  uint64_t bits = 0;
  std::memcpy(&bits, &d, sizeof(double));

  if ((bits >> 63) != 0) {
    return ~bits;
  }

  return bits | 0x8000000000000000ULL;
}

double SortableIntToDouble(uint64_t bits) {
  if ((bits >> 63) != 0) {
    bits &= ~0x8000000000000000ULL;
  } else {
    bits = ~bits;
  }

  double result = 0.0;
  std::memcpy(&result, &bits, sizeof(double));

  return result;
}

size_t NextPowerOfTwo(size_t n) {
  size_t power = 1;
  while (power < n) {
    power <<= 1;
  }
  return power;
}

void RadixSortDouble(std::vector<double> &data) {
  if (data.empty()) {
    return;
  }

  std::vector<uint64_t> keys(data.size());

  for (size_t i = 0; i < data.size(); ++i) {
    keys[i] = DoubleToSortableInt(data[i]);
  }

  std::vector<uint64_t> temp(keys.size());

  constexpr int kRadix = 256;

  for (int pass = 0; pass < 8; ++pass) {
    std::vector<int> count(kRadix, 0);

    int shift = pass * 8;

    for (uint64_t key : keys) {
      count[(key >> shift) & 0xFF]++;
    }

    for (int i = 1; i < kRadix; ++i) {
      count[i] += count[i - 1];
    }

    for (int i = static_cast<int>(keys.size()) - 1; i >= 0; --i) {
      uint8_t byte = (keys[i] >> shift) & 0xFF;
      temp[--count[byte]] = keys[i];
    }

    std::swap(keys, temp);
  }

  for (size_t i = 0; i < data.size(); ++i) {
    data[i] = SortableIntToDouble(keys[i]);
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

void SortChunkTBB(double *raw_data, int chunk_idx, size_t chunk_size) {
  size_t start_idx = static_cast<size_t>(chunk_idx) * chunk_size;

  std::vector<double> local_arr(raw_data + start_idx, raw_data + start_idx + chunk_size);

  RadixSortDouble(local_arr);

  std::ranges::copy(local_arr.begin(), local_arr.end(), raw_data + start_idx);
}

}  // namespace

GonozovLBitSortBatcherMergeALL::GonozovLBitSortBatcherMergeALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());

  GetInput() = in;
}

bool GonozovLBitSortBatcherMergeALL::ValidationImpl() {
  return true;
}

bool GonozovLBitSortBatcherMergeALL::PreProcessingImpl() {
  local_data_ = GetInput();

  return true;
}

bool GonozovLBitSortBatcherMergeALL::RunImpl() {
  int mpi_initialized = 0;
  MPI_Initialized(&mpi_initialized);

  int rank = 0, size = 1;
  if (mpi_initialized) {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
  }

  bool use_mpi = (mpi_initialized && size > 1);
  std::vector<double> local_data;

  if (use_mpi) {
    // ==================== MPI + TBB ВЕРСИЯ ====================

    // 1. Распределение данных
    size_t global_n = local_data_.size();
    MPI_Bcast(&global_n, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    if (global_n == 0) {
      GetOutput() = std::vector<double>();
      return true;
    }

    // Вычисляем размер порции для каждого процесса
    size_t base_count = global_n / size;
    size_t remainder = global_n % size;

    std::vector<int> send_counts(size), send_displs(size);
    size_t offset = 0;
    for (int i = 0; i < size; ++i) {
      send_counts[i] = static_cast<int>(base_count + (static_cast<size_t>(i) < remainder ? 1 : 0));
      send_displs[i] = static_cast<int>(offset);
      offset += send_counts[i];
    }

    int local_n = send_counts[rank];
    local_data.resize(local_n);

    // Рассылка данных
    if (rank == 0) {
      MPI_Scatterv(local_data_.data(), send_counts.data(), send_displs.data(), MPI_DOUBLE, local_data.data(), local_n,
                   MPI_DOUBLE, 0, MPI_COMM_WORLD);
      local_data_.clear();
      local_data_.shrink_to_fit();
    } else {
      MPI_Scatterv(nullptr, nullptr, nullptr, MPI_DOUBLE, local_data.data(), local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    // 2. Локальная TBB-сортировка
    if (!local_data.empty()) {
      size_t local_n_original = local_data.size();
      size_t local_new_size = NextPowerOfTwo(local_n_original);

      if (local_new_size > local_n_original) {
        local_data.resize(local_new_size, std::numeric_limits<double>::infinity());
      }

      int num_threads = ppc::util::GetNumThreads();
      if (num_threads <= 0) {
        num_threads = 1;
      }

      size_t num_chunks = 1;
      while (num_chunks * 2 <= static_cast<size_t>(num_threads) && num_chunks * 2 <= local_new_size) {
        num_chunks *= 2;
      }

      size_t chunk_size = std::max<size_t>(1, local_new_size / num_chunks);
      double *raw_data = local_data.data();
      int num_chunks_int = static_cast<int>(num_chunks);

      tbb::parallel_for(0, num_chunks_int, [&](int i) { SortChunkTBB(raw_data, i, chunk_size); });

      for (size_t cur_size = chunk_size; cur_size < local_new_size; cur_size *= 2) {
        int merges_count = static_cast<int>(local_new_size / (cur_size * 2));
        tbb::parallel_for(0, merges_count, [&](int i) {
          OddEvenMergeIterative(raw_data, static_cast<size_t>(i) * 2 * cur_size, 2 * cur_size);
        });
      }

      if (local_new_size > local_n_original) {
        local_data.resize(local_n_original);
      }
    }

    // 3. Сбор всех данных на процесс 0 (упрощённое слияние)
    if (rank == 0) {
      std::vector<double> all_data = std::move(local_data);
      for (int p = 1; p < size; ++p) {
        int recv_size;
        MPI_Recv(&recv_size, 1, MPI_INT, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (recv_size > 0) {
          std::vector<double> recv_data(recv_size);
          MPI_Recv(recv_data.data(), recv_size, MPI_DOUBLE, p, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          all_data.insert(all_data.end(), recv_data.begin(), recv_data.end());
        }
      }
      local_data = std::move(all_data);
    } else {
      int my_size = static_cast<int>(local_data.size());
      MPI_Send(&my_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      if (my_size > 0) {
        MPI_Send(local_data.data(), my_size, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
      }
      local_data.clear();
    }

    // 4. Финальная сортировка на процессе 0
    if (rank == 0 && !local_data.empty()) {
      RadixSortDouble(local_data);
    }

    // 5. Рассылка результата всем процессам
    size_t total_size = local_data.size();
    MPI_Bcast(&total_size, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);

    if (rank == 0) {
      if (total_size > 0) {
        MPI_Bcast(local_data.data(), static_cast<int>(total_size), MPI_DOUBLE, 0, MPI_COMM_WORLD);
      }
      GetOutput() = std::move(local_data);
    } else {
      local_data.resize(total_size);
      if (total_size > 0) {
        MPI_Bcast(local_data.data(), static_cast<int>(total_size), MPI_DOUBLE, 0, MPI_COMM_WORLD);
      }
      GetOutput() = std::move(local_data);
    }

  } else {
    // ==================== TBB ВЕРСИЯ (БЕЗ MPI) ====================
    local_data = std::move(local_data_);

    if (!local_data.empty()) {
      size_t n = local_data.size();
      size_t new_size = NextPowerOfTwo(n);

      if (new_size > n) {
        local_data.resize(new_size, std::numeric_limits<double>::infinity());
      }

      int num_threads = ppc::util::GetNumThreads();
      if (num_threads <= 0) {
        num_threads = 1;
      }

      size_t num_chunks = 1;
      while (num_chunks * 2 <= static_cast<size_t>(num_threads) && num_chunks * 2 <= new_size) {
        num_chunks *= 2;
      }

      size_t chunk_size = std::max<size_t>(1, new_size / num_chunks);
      double *raw_data = local_data.data();
      int num_chunks_int = static_cast<int>(num_chunks);

      tbb::parallel_for(0, num_chunks_int, [&](int i) { SortChunkTBB(raw_data, i, chunk_size); });

      for (size_t size_step = chunk_size; size_step < new_size; size_step *= 2) {
        int merges_count = static_cast<int>(new_size / (size_step * 2));
        tbb::parallel_for(0, merges_count, [&](int i) {
          OddEvenMergeIterative(raw_data, static_cast<size_t>(i) * 2 * size_step, 2 * size_step);
        });
      }

      if (new_size > n) {
        local_data.resize(n);
      }
    }

    GetOutput() = std::move(local_data);
  }

  return true;
}

bool GonozovLBitSortBatcherMergeALL::PostProcessingImpl() {
  return true;
}

}  // namespace gonozov_l_bitwise_sorting_double_batcher_merge

// #include "gonozov_l_bitwise_sorting_double_Batcher_merge/all/include/ops_all.hpp"

// #include <mpi.h>
// #include <tbb/tbb.h>

// #include <algorithm>
// #include <cstdint>
// #include <cstring>
// #include <limits>
// #include <vector>

// #include "gonozov_l_bitwise_sorting_double_Batcher_merge/common/include/common.hpp"
// #include "util/include/util.hpp"

// namespace gonozov_l_bitwise_sorting_double_batcher_merge {

// namespace {

// uint64_t DoubleToSortableInt(double d) {
//   uint64_t bits = 0;
//   std::memcpy(&bits, &d, sizeof(double));

//   if ((bits >> 63) != 0) {
//     return ~bits;
//   }

//   return bits | 0x8000000000000000ULL;
// }

// double SortableIntToDouble(uint64_t bits) {
//   if ((bits >> 63) != 0) {
//     bits &= ~0x8000000000000000ULL;
//   } else {
//     bits = ~bits;
//   }

//   double result = 0.0;
//   std::memcpy(&result, &bits, sizeof(double));

//   return result;
// }

// size_t NextPowerOfTwo(size_t n) {
//   size_t power = 1;
//   while (power < n) {
//     power <<= 1;
//   }
//   return power;
// }

// void RadixSortDouble(std::vector<double> &data) {
//   if (data.empty()) {
//     return;
//   }

//   std::vector<uint64_t> keys(data.size());

//   for (size_t i = 0; i < data.size(); ++i) {
//     keys[i] = DoubleToSortableInt(data[i]);
//   }

//   std::vector<uint64_t> temp(keys.size());

//   constexpr int kRadix = 256;

//   for (int pass = 0; pass < 8; ++pass) {
//     std::vector<int> count(kRadix, 0);

//     int shift = pass * 8;

//     for (uint64_t key : keys) {
//       count[(key >> shift) & 0xFF]++;
//     }

//     for (int i = 1; i < kRadix; ++i) {
//       count[i] += count[i - 1];
//     }

//     for (int i = static_cast<int>(keys.size()) - 1; i >= 0; --i) {
//       uint8_t byte = (keys[i] >> shift) & 0xFF;
//       temp[--count[byte]] = keys[i];
//     }

//     std::swap(keys, temp);
//   }

//   for (size_t i = 0; i < data.size(); ++i) {
//     data[i] = SortableIntToDouble(keys[i]);
//   }
// }

// void CompareExchangeBlocks(double *arr, size_t i, size_t step) {
//   for (size_t k = 0; k < step; ++k) {
//     if (arr[i + k] > arr[i + k + step]) {
//       std::swap(arr[i + k], arr[i + k + step]);
//     }
//   }
// }

// void OddEvenMergeIterative(double *arr, size_t start, size_t n) {
//   if (n <= 1) {
//     return;
//   }

//   size_t step = n / 2;

//   CompareExchangeBlocks(arr, start, step);

//   step /= 2;

//   for (; step > 0; step /= 2) {
//     for (size_t i = step; i < n - step; i += step * 2) {
//       CompareExchangeBlocks(arr, start + i, step);
//     }
//   }
// }

// void SortChunkALL(double *raw_data, int chunk_idx, size_t chunk_size) {
//   size_t start_idx = static_cast<size_t>(chunk_idx) * chunk_size;

//   std::vector<double> local_arr(raw_data + start_idx, raw_data + start_idx + chunk_size);

//   RadixSortDouble(local_arr);

//   std::ranges::copy(local_arr.begin(), local_arr.end(), raw_data + start_idx);
// }

// }  // namespace

// GonozovLBitSortBatcherMergeALL::GonozovLBitSortBatcherMergeALL(const InType &in) {
//   SetTypeOfTask(GetStaticTypeOfTask());

//   GetInput() = in;
// }

// bool GonozovLBitSortBatcherMergeALL::ValidationImpl() {
//   return true;
// }

// bool GonozovLBitSortBatcherMergeALL::PreProcessingImpl() {
//   local_data_ = GetInput();

//   return true;
// }

// bool GonozovLBitSortBatcherMergeALL::RunImpl() {
//   if (local_data_.empty()) {
//     return true;
//   }

//   size_t n = local_data_.size();

//   size_t new_size = NextPowerOfTwo(n);

//   if (new_size > n) {
//     local_data_.resize(new_size, std::numeric_limits<double>::infinity());
//   }

//   int num_threads = ppc::util::GetNumThreads();

//   if (num_threads <= 0) {
//     num_threads = 1;
//   }

//   size_t num_chunks = 1;

//   while (num_chunks * 2 <= static_cast<size_t>(num_threads) && num_chunks * 2 <= new_size) {
//     num_chunks *= 2;
//   }

//   size_t chunk_size = new_size / num_chunks;
//   chunk_size = std::max<size_t>(1, chunk_size);
//   double *raw_data = local_data_.data();
//   int num_chunks_int = static_cast<int>(num_chunks);

//   tbb::parallel_for(0, num_chunks_int, [&](int i) { SortChunkALL(raw_data, i, chunk_size); });

//   for (size_t size = chunk_size; size < new_size; size *= 2) {
//     int merges_count = static_cast<int>(new_size / (size * 2));

//     tbb::parallel_for(0, merges_count,
//                       [&](int i) { OddEvenMergeIterative(raw_data, static_cast<size_t>(i) * 2 * size, 2 * size); });
//   }

//   if (new_size > n) {
//     local_data_.resize(n);
//   }

//   MPI_Barrier(MPI_COMM_WORLD);
//   GetOutput() = local_data_;
//   return true;
// }

// bool GonozovLBitSortBatcherMergeALL::PostProcessingImpl() {
//   return true;
// }

// }  // namespace gonozov_l_bitwise_sorting_double_batcher_merge
