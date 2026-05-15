#include "gonozov_l_bitwise_sorting_double_Batcher_merge/stl/include/ops_stl.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <thread>
#include <vector>

#include "gonozov_l_bitwise_sorting_double_Batcher_merge/common/include/common.hpp"

namespace gonozov_l_bitwise_sorting_double_batcher_merge {

GonozovLBitSortBatcherMergeSTL::GonozovLBitSortBatcherMergeSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool GonozovLBitSortBatcherMergeSTL::ValidationImpl() {
  return !GetInput().empty();  // проверка на то, что исходный массив непустой
}

bool GonozovLBitSortBatcherMergeSTL::PreProcessingImpl() {
  return true;
}

namespace {

constexpr size_t kRadix = 256;
constexpr size_t kCutoff = 1 << 18;

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

void RadixSortDouble(std::vector<double> &data) {
  if (data.empty()) {
    return;
  }

  std::vector<uint64_t> keys(data.size());

  for (size_t i = 0; i < data.size(); ++i) {
    keys[i] = DoubleToSortableInt(data[i]);
  }

  constexpr size_t kByteRange = 256;
  std::vector<uint64_t> temp_keys(data.size());

  for (size_t byte_id = 0; byte_id < 8; ++byte_id) {
    std::array<size_t, kByteRange> count{};

    size_t shift = byte_id * 8;

    for (uint64_t key : keys) {
      uint8_t byte = static_cast<uint8_t>((key >> shift) & 0xFF);
      ++count[byte];
    }

    for (size_t i = 1; i < kByteRange; ++i) {
      count[i] += count[i - 1];
    }

    for (size_t i = keys.size(); i-- > 0;) {
      uint8_t byte = static_cast<uint8_t>((keys[i] >> shift) & 0xFF);
      temp_keys[--count[byte]] = keys[i];
    }

    keys.swap(temp_keys);
  }

  for (size_t i = 0; i < data.size(); ++i) {
    data[i] = SortableIntToDouble(keys[i]);
  }
}

void SortChunk(double *data, size_t start, size_t size) {
  std::vector<double> local(data + start, data + start + size);

  RadixSortDouble(local);

  std::ranges::copy(local, data + static_cast<ptrdiff_t>(start));
}

void MergeSortedChunks(double *data, size_t total_size, size_t chunk_size) {
  for (size_t current_size = chunk_size; current_size < total_size; current_size *= 2) {
    for (size_t left = 0; left < total_size; left += current_size * 2) {
      std::inplace_merge(
          data + left,
          data + left + current_size,
          data + left + current_size * 2);
    }
  }
}

size_t NextPowerOfTwo(size_t n) {
  size_t power = 1;

  while (power < n) {
    power <<= 1;
  }

  return power;
}

void ParallelHybridSort(std::vector<double> &data) {
  if (data.size() <= 1) {
    return;
  }

  size_t original_size = data.size();
  size_t padded_size = NextPowerOfTwo(original_size);

  data.resize(padded_size, std::numeric_limits<double>::infinity());

  int threads_count = ppc::util::GetNumThreads();

  threads_count = std::max(1, threads_count);

  size_t chunks_count = 1;

  while (chunks_count * 2 <= static_cast<size_t>(threads_count) &&
         chunks_count * 2 <= padded_size) {
    chunks_count *= 2;
  }

  size_t chunk_size = padded_size / chunks_count;

  double *raw_data = data.data();

  std::vector<std::thread> threads;
  threads.reserve(chunks_count);

  for (size_t i = 0; i < chunks_count; ++i) {
    threads.emplace_back([raw_data, i, chunk_size]() {
      SortChunk(raw_data, i * chunk_size, chunk_size);
    });
  }

  for (auto &thread : threads) {
    thread.join();
  }

  MergeSortedChunks(raw_data, padded_size, chunk_size);

  data.resize(original_size);
}

// uint64_t DoubleToUint(double d) {
//   uint64_t u = 0;
//   std::memcpy(&u, &d, sizeof(double));
//   if ((u & 0x8000000000000000ULL) != 0) {
//     u = ~u;
//   } else {
//     u |= 0x8000000000000000ULL;
//   }
//   return u;
// }

// double UintToDouble(uint64_t u) {
//   if ((u & 0x8000000000000000ULL) != 0) {
//     u &= ~0x8000000000000000ULL;
//   } else {
//     u = ~u;
//   }
//   double d = 0.0;
//   std::memcpy(&d, &u, sizeof(double));
//   return d;
// }

// void RadixSortDouble(std::vector<double> &arr) {
//   if (arr.empty()) {
//     return;
//   }

//   std::vector<uint64_t> uarr(arr.size());
//   for (size_t i = 0; i < arr.size(); ++i) {
//     uarr[i] = DoubleToUint(arr[i]);
//   }

//   std::vector<uint64_t> temp(uarr.size());
//   for (size_t byte = 0; byte < 8; ++byte) {
//   std::array<size_t, 256> count{};

//   for (uint64_t val : uarr) {
//     count[(val >> (byte * 8)) & 0xFF]++;
//   }

//   for (size_t i = 1; i < 256; ++i) {
//     count[i] += count[i - 1];
//   }

//   for (size_t i = uarr.size(); i-- > 0;) {
//     temp[--count[(uarr[i] >> (byte * 8)) & 0xFF]] = uarr[i];
//   }

//   uarr = temp;
// }

//   for (size_t i = 0; i < arr.size(); ++i) {
//     arr[i] = UintToDouble(uarr[i]);
//   }
// }

// void ProcessChunkSTL(double *raw_data, int chunk_idx, size_t chunk_size) {
//   size_t start_idx = static_cast<size_t>(chunk_idx) * chunk_size;

//   std::vector<double> local_arr(chunk_size);
//   double *local_raw = local_arr.data();

//   for (size_t j = 0; j < chunk_size; ++j) {
//     local_raw[j] = raw_data[start_idx + j];
//   }

//   RadixSortDouble(local_arr);

//   for (size_t j = 0; j < chunk_size; ++j) {
//     raw_data[start_idx + j] = local_raw[j];
//   }
// }

// void ExecuteSTLSort(double *raw_data, size_t pow2, size_t chunk_size, int num_chunks_int) {
//   std::vector<std::thread> threads;
//   threads.reserve(num_chunks_int);

//   // 1. Параллельно сортируем каждый блок, раздавая задачи сырым std::thread
//   for (int i = 0; i < num_chunks_int; ++i) {
//     threads.emplace_back([raw_data, i, chunk_size]() { ProcessChunkSTL(raw_data, i, chunk_size); });
//   }

//   for (auto &t : threads) {
//     t.join();
//   }
//   threads.clear();

// // 2. Параллельное слияние отсортированных блоков
// for (size_t size = chunk_size; size < pow2; size *= 2) {
//   for (size_t start = 0; start < pow2; start += 2 * size) {
//     std::inplace_merge(
//         raw_data + start,
//         raw_data + start + size,
//         raw_data + start + 2 * size);
//   }
// }
// }

}// namespace

bool GonozovLBitSortBatcherMergeSTL::RunImpl() {
  std::vector<double> array = GetInput();

  ParallelHybridSort(array);

  GetOutput() = array;
  // size_t original_size = array.size();
  // size_t pow2 = 1;
  // while (pow2 < original_size) {
  //   pow2 *= 2;
  // }

  // if (pow2 > original_size) {
  //   array.resize(pow2, std::numeric_limits<double>::infinity());
  // }

  // int num_threads = ppc::util::GetNumThreads();
  // if (num_threads <= 0) {
  //   num_threads = 1;
  // }

  // size_t num_chunks = 1;
  // while (num_chunks * 2 <= static_cast<size_t>(num_threads) && num_chunks * 2 <= pow2) {
  //   num_chunks *= 2;
  // }

  // size_t chunk_size = pow2 / num_chunks;
  // double *raw_data = array.data();
  // int num_chunks_int = static_cast<int>(num_chunks);

  // ExecuteSTLSort(raw_data, pow2, chunk_size, num_chunks_int);

  // if (pow2 > original_size) {
  //   array.resize(original_size);
  // }

  // GetOutput() = array;
  return true;
}

bool GonozovLBitSortBatcherMergeSTL::PostProcessingImpl() {
  return true;
}

}  // namespace gonozov_l_bitwise_sorting_double_batcher_merge
