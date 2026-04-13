#include "papulina_y_radix_sort/omp/include/ops_omp.hpp"

#include <omp.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <span>
#include <utility>
#include <vector>

#include "papulina_y_radix_sort/common/include/common.hpp"

namespace papulina_y_radix_sort {

PapulinaYRadixSortOMP::PapulinaYRadixSortOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<double>();
}
bool PapulinaYRadixSortOMP::ValidationImpl() {
  return true;
}

bool PapulinaYRadixSortOMP::PreProcessingImpl() {
  return true;
}
bool PapulinaYRadixSortOMP::RunImpl() {
  double *result = GetInput().data();
  int size = static_cast<int>(GetInput().size());
  // int threads_count = std::min(omp_get_max_threads(), std::max(4, size / 1000));
  int threads_count = omp_get_max_threads();

  std::vector<std::span<double>> chunks;
  std::vector<int> chunks_offsets;

  int chunk_size = size / threads_count;
  int reminder = size % threads_count;

  int offset = 0;
  for (int i = 0; i < threads_count; i++) {
    int real_chunk_size = chunk_size + (i < reminder ? 1 : 0);
    if (real_chunk_size > 0) {
      chunks_offsets.push_back(offset);
      chunks.emplace_back(result + offset, real_chunk_size);
      offset += real_chunk_size;
    }
  }
  threads_count = static_cast<int>(
      chunks.size());  // тк возможно chunks.size() < threads_count(каким-то потокам ничего не распределилось из данных)

#pragma omp parallel for default(none) shared(result, chunks, threads_count) num_threads(threads_count)
  for (int i = 0; i < threads_count; i++) {
    RadixSort(chunks[i].data(), static_cast<int>(chunks[i].size()));
  }

  MergeChunks(chunks, result);
  GetOutput() = std::vector<double>(size);
  for (int i = 0; i < size; i++) {
    GetOutput()[i] = result[i];
    // std::cout << result [i] << " ";
  }
  // std::cout << std::endl;
  return true;
}
std::vector<double> PapulinaYRadixSortOMP::Merge(const std::vector<double> &a, const std::vector<double> &b) {
  std::vector<double> result;
  result.reserve(a.size() + b.size());
  size_t i = 0;
  size_t j = 0;
  while (i < a.size() && j < b.size()) {
    if (a[i] <= b[j]) {
      result.push_back(a[i]);
      ++i;
    } else {
      result.push_back(b[j]);
      ++j;
    }
  }
  while (i < a.size()) {
    result.push_back(a[i]);
    ++i;
  }
  while (j < b.size()) {
    result.push_back(b[j]);
    ++j;
  }
  return result;
}
void PapulinaYRadixSortOMP::MergeChunks(std::vector<std::span<double>> &chunks, double *result) {
  if (chunks.size() <= 1) {
    return;
  }

  int n = static_cast<int>(GetInput().size());

  std::vector<std::vector<double>> chunks_copy;
  chunks_copy.reserve(chunks.size());
  for (auto &chunk : chunks) {
    chunks_copy.emplace_back(chunk.begin(), chunk.end());
  }
  while (chunks_copy.size() > 1) {
    size_t pair_count = chunks_copy.size() / 2;
    std::vector<std::vector<double>> next((chunks_copy.size() + 1) / 2);

#pragma omp parallel for default(none) shared(chunks_copy, next, pair_count)
    for (size_t i = 0; i < pair_count; ++i) {
      size_t idx = i;
      next[idx] = Merge(chunks_copy[2 * idx], chunks_copy[(2 * idx) + 1]);
    }
    if (chunks_copy.size() % 2 != 0) {
      next.back() = std::move(chunks_copy.back());
    }
    chunks_copy = std::move(next);
  }
  const auto &final_result = chunks_copy[0];
  for (int i = 0; i < n; ++i) {
    result[i] = final_result[i];
  }
}
void PapulinaYRadixSortOMP::RadixSort(double *arr, int size) {
  std::vector<uint64_t> bytes(size);
  std::vector<uint64_t> out(size);

  for (int i = 0; i < size; i++) {
    bytes[i] = InBytes(arr[i]);
  }

  SortByByte(bytes.data(), out.data(), 0, size);
  SortByByte(out.data(), bytes.data(), 1, size);
  SortByByte(bytes.data(), out.data(), 2, size);
  SortByByte(out.data(), bytes.data(), 3, size);
  SortByByte(bytes.data(), out.data(), 4, size);
  SortByByte(out.data(), bytes.data(), 5, size);
  SortByByte(bytes.data(), out.data(), 6, size);
  SortByByte(out.data(), bytes.data(), 7, size);

  for (int i = 0; i < size; i++) {
    arr[i] = FromBytes(bytes[i]);
  }
}
bool PapulinaYRadixSortOMP::PostProcessingImpl() {
  return true;
}
void PapulinaYRadixSortOMP::SortByByte(uint64_t *bytes, uint64_t *out, int byte, int size) {
  auto *byte_view = reinterpret_cast<unsigned char *>(bytes);  // просматриваем как массив байтов
  std::array<int, 256> counter = {0};

  for (int i = 0; i < size; i++) {
    int index = byte_view[(8 * i) + byte];
    *(counter.data() + index) += 1;
  }
  int tmp = 0;
  int j = 0;
  for (; j < 256; j++) {
    if (*(counter.data() + j) != 0) {
      tmp = *(counter.data() + j);
      *(counter.data() + j) = 0;
      j++;
      break;
    }
  }
  for (; j < 256; j++) {
    int a = *(counter.data() + j);
    *(counter.data() + j) = tmp;
    tmp += a;
  }
  for (int i = 0; i < size; i++) {
    int index = byte_view[(8 * i) + byte];
    out[*(counter.data() + index)] = bytes[i];
    *(counter.data() + index) += 1;
  }
}
uint64_t PapulinaYRadixSortOMP::InBytes(double d) {
  uint64_t bits = 0;
  memcpy(&bits, &d, sizeof(double));
  if ((bits & kMask) != 0) {
    bits = ~bits;
  } else {
    bits = bits ^ kMask;
  }
  return bits;
}
double PapulinaYRadixSortOMP::FromBytes(uint64_t bits) {
  double d = NAN;
  if ((bits & kMask) != 0) {
    bits = bits ^ kMask;
  } else {
    bits = ~bits;
  }
  memcpy(&d, &bits, sizeof(double));
  return d;
}
}  // namespace papulina_y_radix_sort
