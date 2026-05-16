#include "frolova_s_radix_sort_double/all/include/ops_all.hpp"

#include <omp.h>
#include <tbb/parallel_for.h>
#include <tbb/task_group.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>

namespace frolova_s_radix_sort_double {

FrolovaSRadixSortDoubleALL::FrolovaSRadixSortDoubleALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

uint64_t FrolovaSRadixSortDoubleALL::InBytes(double d) {
  uint64_t bits;
  std::memcpy(&bits, &d, sizeof(double));
  if ((bits & kMask) != 0) {
    bits = ~bits;
  } else {
    bits = bits ^ kMask;
  }
  return bits;
}

double FrolovaSRadixSortDoubleALL::FromBytes(uint64_t bits) {
  if ((bits & kMask) != 0) {
    bits = bits ^ kMask;
  } else {
    bits = ~bits;
  }
  double d;
  std::memcpy(&d, &bits, sizeof(double));
  return d;
}

void FrolovaSRadixSortDoubleALL::SortByByte(uint64_t *bytes, uint64_t *out, int byte, int size) {
  auto *byte_view = reinterpret_cast<unsigned char *>(bytes);
  std::array<int, 256> counter = {0};

  for (int i = 0; i < size; i++) {
    int index = byte_view[(8 * i) + byte];
    counter[index]++;
  }

  int total = 0;
  for (int j = 0; j < 256; j++) {
    int old = counter[j];
    counter[j] = total;
    total += old;
  }

  for (int i = 0; i < size; i++) {
    int index = byte_view[(8 * i) + byte];
    out[counter[index]] = bytes[i];
    counter[index]++;
  }
}

void FrolovaSRadixSortDoubleALL::RadixSort(double *arr, int size) {
  if (size <= 1) {
    return;
  }

  std::vector<uint64_t> bytes(size);
  std::vector<uint64_t> out(size);

#pragma omp parallel for
  for (int i = 0; i < size; i++) {
    bytes[i] = InBytes(arr[i]);
  }

  uint64_t *src_ptr = bytes.data();
  uint64_t *dst_ptr = out.data();

  for (int byte = 0; byte < 8; byte++) {
    SortByByte(src_ptr, dst_ptr, byte, size);
    std::swap(src_ptr, dst_ptr);
  }

#pragma omp parallel for
  for (int i = 0; i < size; i++) {
    arr[i] = FromBytes(src_ptr[i]);
  }
}

std::vector<double> FrolovaSRadixSortDoubleALL::SimpleMerge(const std::vector<double> &a,
                                                            const std::vector<double> &b) {
  std::vector<double> res;
  res.reserve(a.size() + b.size());
  size_t i = 0, j = 0;
  while (i < a.size() && j < b.size()) {
    if (a[i] <= b[j]) {
      res.push_back(a[i++]);
    } else {
      res.push_back(b[j++]);
    }
  }
  while (i < a.size()) {
    res.push_back(a[i++]);
  }
  while (j < b.size()) {
    res.push_back(b[j++]);
  }
  return res;
}

std::vector<double> FrolovaSRadixSortDoubleALL::ParallelMerge(std::vector<std::vector<double>> &chunks) {
  if (chunks.empty()) {
    return {};
  }
  if (chunks.size() == 1) {
    return std::move(chunks[0]);
  }

  tbb::task_group tg;
  std::vector<std::vector<double>> next_chunks;
  next_chunks.reserve((chunks.size() + 1) / 2);

  for (size_t i = 0; i < chunks.size(); i += 2) {
    if (i + 1 < chunks.size()) {
      next_chunks.emplace_back();
      tg.run([&, i] { next_chunks.back() = SimpleMerge(chunks[i], chunks[i + 1]); });
    } else {
      next_chunks.push_back(std::move(chunks[i]));
    }
  }
  tg.wait();

  return ParallelMerge(next_chunks);
}

bool FrolovaSRadixSortDoubleALL::ValidationImpl() {
  return true;
}
bool FrolovaSRadixSortDoubleALL::PreProcessingImpl() {
  return true;
}

bool FrolovaSRadixSortDoubleALL::RunImpl() {
  std::vector<double> &input = GetInput();
  if (input.empty()) {
    return true;
  }

  int num_chunks = omp_get_max_threads();
  int total_size = static_cast<int>(input.size());
  int chunk_size = total_size / num_chunks;
  int remainder = total_size % num_chunks;

  std::vector<std::vector<double>> chunks(num_chunks);
  int offset = 0;
  for (int i = 0; i < num_chunks; ++i) {
    int cur_size = chunk_size + (i < remainder ? 1 : 0);
    if (cur_size > 0) {
      chunks[i].assign(input.begin() + offset, input.begin() + offset + cur_size);
      offset += cur_size;
    }
  }

  tbb::parallel_for(0, num_chunks, [&](int i) {
    if (!chunks[i].empty()) {
      RadixSort(chunks[i].data(), static_cast<int>(chunks[i].size()));
    }
  });

  std::vector<double> sorted = ParallelMerge(chunks);
  GetOutput() = std::move(sorted);

  return true;
}

bool FrolovaSRadixSortDoubleALL::PostProcessingImpl() {
  return true;
}

}  // namespace frolova_s_radix_sort_double
