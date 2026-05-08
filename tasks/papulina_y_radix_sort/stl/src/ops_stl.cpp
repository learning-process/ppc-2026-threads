#include "papulina_y_radix_sort/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <thread>
#include <utility>
#include <vector>

#include "papulina_y_radix_sort/common/include/common.hpp"
namespace papulina_y_radix_sort {

PapulinaYRadixSortSTL::PapulinaYRadixSortSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<double>();
}

bool PapulinaYRadixSortSTL::ValidationImpl() {
  return !GetInput().empty();
}

bool PapulinaYRadixSortSTL::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return true;
}

bool PapulinaYRadixSortSTL::RunImpl() {
  size_t size = GetInput().size();
  if (size == 0) {
    return true;
  }

  std::vector<double> data = GetInput();
  RadixSortParallel(data.data(), size);

  GetOutput() = std::move(data);
  return true;
}

bool PapulinaYRadixSortSTL::PostProcessingImpl() {
  return true;
}

uint64_t PapulinaYRadixSortSTL::InBytes(double d) {
  uint64_t bits = 0;
  memcpy(&bits, &d, sizeof(double));
  if ((bits & kMask) != 0) {
    bits = ~bits;
  } else {
    bits = bits ^ kMask;
  }
  return bits;
}

double PapulinaYRadixSortSTL::FromBytes(uint64_t bits) {
  double d = 0;
  if ((bits & kMask) != 0) {
    bits = bits ^ kMask;
  } else {
    bits = ~bits;
  }
  memcpy(&d, &bits, sizeof(double));
  return d;
}

void PapulinaYRadixSortSTL::RadixSortParallel(double *arr, size_t size) {
  if (arr == nullptr || size == 0) {
    return;
  }

  unsigned int num_threads = std::max(1U, std::thread::hardware_concurrency());
  if (size < 1000) {
    num_threads = 3;
  }

  std::vector<uint64_t> bytes;
  bytes.resize(size);
  for (size_t i = 0; i < size; ++i) {
    bytes[i] = InBytes(arr[i]);
  }

  std::vector<uint64_t> temp;
  temp.resize(size);

  uint64_t *src_ptr = bytes.data();
  uint64_t *dst_ptr = temp.data();

  for (int byte_idx = 0; byte_idx < 8; ++byte_idx) {
    ExecuteRadixPass(src_ptr, dst_ptr, size, byte_idx, num_threads);
    std::swap(src_ptr, dst_ptr);
  }

  for (size_t i = 0; i < size; ++i) {
    arr[i] = FromBytes(src_ptr[i]);
  }
}

std::pair<size_t, size_t> PapulinaYRadixSortSTL::GetThreadRange(size_t size, unsigned int num_threads,
                                                                unsigned int thread_idx) {
  size_t chunk = size / num_threads;
  size_t start = thread_idx * chunk;
  size_t end = (thread_idx == num_threads - 1) ? size : (thread_idx + 1) * chunk;
  return {start, end};
}

void PapulinaYRadixSortSTL::CountFrequencies(const uint64_t *src, size_t size, int byte_idx, unsigned int num_threads,
                                             std::vector<std::vector<size_t>> &local_hists) {
  std::vector<std::thread> workers;
  workers.reserve(num_threads);
  for (unsigned int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    workers.emplace_back([&, thread_idx]() {
      auto range = GetThreadRange(size, num_threads, thread_idx);
      const auto *byte_view = reinterpret_cast<const unsigned char *>(src);
      for (size_t i = range.first; i < range.second; ++i) {
        local_hists[thread_idx][byte_view[(i * 8) + byte_idx]]++;
      }
    });
  }
  for (auto &worker : workers) {
    worker.join();
  }
}

void PapulinaYRadixSortSTL::ComputeOffsets(unsigned int num_threads,
                                           const std::vector<std::vector<size_t>> &local_hists,
                                           std::vector<std::vector<size_t>> &thread_pos) {
  size_t total = 0;
  for (int bucket = 0; bucket < 256; ++bucket) {
    for (unsigned int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
      thread_pos[thread_idx][bucket] = total;
      total += local_hists[thread_idx][bucket];
    }
  }
}

void PapulinaYRadixSortSTL::ReorderElements(const uint64_t *src, uint64_t *dst, size_t size, int byte_idx,
                                            unsigned int num_threads, std::vector<std::vector<size_t>> &thread_pos) {
  std::vector<std::thread> workers;
  workers.reserve(num_threads);
  for (unsigned int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    workers.emplace_back([&, thread_idx]() {
      auto range = GetThreadRange(size, num_threads, thread_idx);
      const auto *byte_view = reinterpret_cast<const unsigned char *>(src);
      for (size_t i = range.first; i < range.second; ++i) {
        int bucket = byte_view[(i * 8) + byte_idx];
        dst[thread_pos[thread_idx][bucket]++] = src[i];
      }
    });
  }
  for (auto &worker : workers) {
    worker.join();
  }
}

void PapulinaYRadixSortSTL::ExecuteRadixPass(uint64_t *src, uint64_t *dst, size_t size, int byte_idx,
                                             unsigned int num_threads) {
  std::vector<std::vector<size_t>> local_hists(num_threads, std::vector<size_t>(256, 0));
  std::vector<std::vector<size_t>> thread_pos(num_threads, std::vector<size_t>(256, 0));

  CountFrequencies(src, size, byte_idx, num_threads, local_hists);
  ComputeOffsets(num_threads, local_hists, thread_pos);
  ReorderElements(src, dst, size, byte_idx, num_threads, thread_pos);
}

}  // namespace papulina_y_radix_sort
