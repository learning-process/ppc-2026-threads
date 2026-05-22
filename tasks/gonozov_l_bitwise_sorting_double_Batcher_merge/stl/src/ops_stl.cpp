// #include "gonozov_l_bitwise_sorting_double_Batcher_merge/stl/include/ops_stl.hpp"

// #include <algorithm>
// #include <cmath>
// #include <cstddef>
// #include <cstdint>
// #include <cstring>
// #include <limits>
// #include <thread>
// #include <vector>

// #include "gonozov_l_bitwise_sorting_double_Batcher_merge/common/include/common.hpp"
// #include "util/include/util.hpp"

// namespace gonozov_l_bitwise_sorting_double_batcher_merge {

// GonozovLBitSortBatcherMergeSTL::GonozovLBitSortBatcherMergeSTL(const InType &in) {
//   SetTypeOfTask(GetStaticTypeOfTask());
//   GetInput() = in;
// }

// bool GonozovLBitSortBatcherMergeSTL::ValidationImpl() {
//   return !GetInput().empty();  // проверка на то, что исходный массив непустой
// }

// bool GonozovLBitSortBatcherMergeSTL::PreProcessingImpl() {
//   return true;
// }

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

// void RadixSortDouble(std::vector<double> &data) {
//   if (data.empty()) {
//     return;
//   }

//   std::vector<uint64_t> keys(data.size());

//   for (size_t i = 0; i < data.size(); ++i) {
//     keys[i] = DoubleToSortableInt(data[i]);
//   }

//   constexpr size_t kByteRange = 256;
//   std::vector<uint64_t> temp_keys(data.size());

//   for (size_t byte_id = 0; byte_id < 8; ++byte_id) {
//     std::vector<size_t> count(kByteRange, 0);

//     size_t shift = byte_id * 8;

//     for (uint64_t key : keys) {
//       auto byte = static_cast<uint8_t>((key >> shift) & 0xFF);
//       ++count[byte];
//     }

//     for (size_t i = 1; i < kByteRange; ++i) {
//       count[i] += count[i - 1];
//     }

//     for (size_t i = keys.size(); i-- > 0;) {
//       auto byte = static_cast<uint8_t>((keys[i] >> shift) & 0xFF);
//       temp_keys[--count[byte]] = keys[i];
//     }

//     keys.swap(temp_keys);
//   }

//   for (size_t i = 0; i < data.size(); ++i) {
//     data[i] = SortableIntToDouble(keys[i]);
//   }
// }

// void SortChunk(double *data, size_t start, size_t size) {
//   std::vector<double> local(data + start, data + start + size);

//   RadixSortDouble(local);

//   std::ranges::copy(local, data + static_cast<ptrdiff_t>(start));
// }

// void MergeSortedChunks(double *data, size_t total_size, size_t chunk_size) {
//   for (size_t current_size = chunk_size; current_size < total_size; current_size *= 2) {
//     for (size_t left = 0; left < total_size; left += current_size * 2) {
//       std::inplace_merge(data + left, data + left + current_size, data + left + (current_size * 2));
//     }
//   }
// }

// size_t NextPowerOfTwo(size_t n) {
//   size_t power = 1;

//   while (power < n) {
//     power <<= 1;
//   }

//   return power;
// }

// void ParallelHybridSort(std::vector<double> &data) {
//   if (data.size() <= 1) {
//     return;
//   }

//   size_t original_size = data.size();
//   size_t padded_size = NextPowerOfTwo(original_size);

//   data.resize(padded_size, std::numeric_limits<double>::infinity());

//   int threads_count = ppc::util::GetNumThreads();

//   threads_count = std::max(1, threads_count);

//   size_t chunks_count = 1;

//   while (chunks_count * 2 <= static_cast<size_t>(threads_count) && chunks_count * 2 <= padded_size) {
//     chunks_count *= 2;
//   }

//   size_t chunk_size = padded_size / chunks_count;

//   double *raw_data = data.data();

//   std::vector<std::thread> threads;
//   threads.reserve(chunks_count);

//   for (size_t i = 0; i < chunks_count; ++i) {
//     threads.emplace_back([raw_data, i, chunk_size]() { SortChunk(raw_data, i * chunk_size, chunk_size); });
//   }

//   for (auto &thread : threads) {
//     thread.join();
//   }

//   MergeSortedChunks(raw_data, padded_size, chunk_size);

//   data.resize(original_size);
// }

// }  // namespace

// bool GonozovLBitSortBatcherMergeSTL::RunImpl() {
//   std::vector<double> array = GetInput();

//   ParallelHybridSort(array);

//   GetOutput() = array;

//   return true;
// }

// bool GonozovLBitSortBatcherMergeSTL::PostProcessingImpl() {
//   return true;
// }

// }  // namespace gonozov_l_bitwise_sorting_double_batcher_merge

#include "gonozov_l_bitwise_sorting_double_Batcher_merge/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <thread>
#include <vector>

#include "gonozov_l_bitwise_sorting_double_Batcher_merge/common/include/common.hpp"
#include "util/include/util.hpp"

namespace gonozov_l_bitwise_sorting_double_batcher_merge {

GonozovLBitSortBatcherMergeSTL::GonozovLBitSortBatcherMergeSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool GonozovLBitSortBatcherMergeSTL::ValidationImpl() {
  return !GetInput().empty();
}

bool GonozovLBitSortBatcherMergeSTL::PreProcessingImpl() {
  return true;
}

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
    std::vector<size_t> count(kByteRange, 0);

    size_t shift = byte_id * 8;

    for (uint64_t key : keys) {
      auto byte = static_cast<uint8_t>((key >> shift) & 0xFF);
      ++count[byte];
    }

    for (size_t i = 1; i < kByteRange; ++i) {
      count[i] += count[i - 1];
    }

    for (size_t i = keys.size(); i-- > 0;) {
      auto byte = static_cast<uint8_t>((keys[i] >> shift) & 0xFF);
      temp_keys[--count[byte]] = keys[i];
    }

    keys.swap(temp_keys);
  }

  for (size_t i = 0; i < data.size(); ++i) {
    data[i] = SortableIntToDouble(keys[i]);
  }
}

// Функция теперь принимает и возвращает копию вектора, а не работает со ссылкой на данные
std::vector<double> SortChunk(const std::vector<double> &chunk) {
  std::vector<double> result = chunk;  // Создаём копию
  RadixSortDouble(result);
  return result;  // Возвращаем отсортированную копию
}

// Функция слияния принимает копии отсортированных блоков и возвращает результат
std::vector<double> MergeTwoChunks(const std::vector<double> &left, const std::vector<double> &right) {
  std::vector<double> result;
  result.reserve(left.size() + right.size());
  result.insert(result.end(), left.begin(), left.end());
  result.insert(result.end(), right.begin(), right.end());

  // Применяем std::inplace_merge, но так как result — это отдельная копия, всё безопасно
  std::inplace_merge(result.begin(), result.begin() + left.size(), result.end());

  return result;
}

// Рекурсивное слияние всех блоков с копированием
std::vector<double> MergeAllChunks(const std::vector<std::vector<double>> &chunks) {
  if (chunks.empty()) {
    return {};
  }

  if (chunks.size() == 1) {
    return chunks[0];
  }

  std::vector<std::vector<double>> merged_chunks;
  merged_chunks.reserve((chunks.size() + 1) / 2);

  for (size_t i = 0; i < chunks.size(); i += 2) {
    if (i + 1 < chunks.size()) {
      merged_chunks.push_back(MergeTwoChunks(chunks[i], chunks[i + 1]));
    } else {
      merged_chunks.push_back(chunks[i]);  // Нечётный блок остаётся без пары
    }
  }

  return MergeAllChunks(merged_chunks);  // Рекурсивно продолжаем слияние
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

  // Создаём расширенный вектор с копиями элементов
  std::vector<double> padded_data = data;
  padded_data.resize(padded_size, std::numeric_limits<double>::infinity());

  int threads_count = ppc::util::GetNumThreads();
  threads_count = std::max(1, threads_count);

  size_t chunks_count = 1;
  while (chunks_count * 2 <= static_cast<size_t>(threads_count) && chunks_count * 2 <= padded_size) {
    chunks_count *= 2;
  }

  size_t chunk_size = padded_size / chunks_count;

  // Разбиваем данные на чанки (копирование элементов)
  std::vector<std::vector<double>> chunks(chunks_count);
  for (size_t i = 0; i < chunks_count; ++i) {
    size_t start = i * chunk_size;
    chunks[i].assign(padded_data.begin() + static_cast<ptrdiff_t>(start),
                     padded_data.begin() + static_cast<ptrdiff_t>(start + chunk_size));
  }

  // Параллельная сортировка чанков (каждый поток работает со своей копией)
  std::vector<std::thread> threads;
  threads.reserve(chunks_count);

  std::vector<std::vector<double>> sorted_chunks(chunks_count);

  for (size_t i = 0; i < chunks_count; ++i) {
    threads.emplace_back([&sorted_chunks, &chunks, i]() {
      sorted_chunks[i] = SortChunk(chunks[i]);  // Каждый поток записывает результат в свою позицию
    });
  }

  for (auto &thread : threads) {
    thread.join();
  }

  // Слияние всех отсортированных чанков с копированием
  std::vector<double> result = MergeAllChunks(sorted_chunks);

  // Копируем результат обратно в исходные данные (обрезая до исходного размера)
  data.assign(result.begin(), result.begin() + static_cast<ptrdiff_t>(original_size));
}

}  // namespace

bool GonozovLBitSortBatcherMergeSTL::RunImpl() {
  std::vector<double> array = GetInput();
  ParallelHybridSort(array);
  GetOutput() = array;
  return true;
}

bool GonozovLBitSortBatcherMergeSTL::PostProcessingImpl() {
  return true;
}

}  // namespace gonozov_l_bitwise_sorting_double_batcher_merge
