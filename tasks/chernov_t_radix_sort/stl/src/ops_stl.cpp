#include "chernov_t_radix_sort/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <thread>
#include <vector>

#include "chernov_t_radix_sort/common/include/common.hpp"
#include "util/include/util.hpp"

namespace chernov_t_radix_sort {

ChernovTRadixSortSTL::ChernovTRadixSortSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool ChernovTRadixSortSTL::ValidationImpl() {
  return true;
}

bool ChernovTRadixSortSTL::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

constexpr int kBitsPerDigit = 8;
constexpr int kRadix = 1 << kBitsPerDigit;  // 256
constexpr uint32_t kSignMask = 0x80000000U;

void ChernovTRadixSortSTL::RadixSortLSDParallel(std::vector<int> &data, int num_threads) {
  const size_t n = data.size();
  if (n <= 1) {
    return;
  }

  std::vector<uint32_t> temp(n);
  std::vector<std::thread> threads;
  const size_t chunk_size = (n + num_threads - 1) / num_threads;

  for (int t = 0; t < num_threads; ++t) {
    const size_t start = t * chunk_size;
    const size_t end = std::min(start + chunk_size, n);
    threads.emplace_back([&data, &temp, start, end]() {
      for (size_t i = start; i < end; ++i) {
        temp[i] = static_cast<uint32_t>(data[i]) ^ kSignMask;
      }
    });
  }
  for (auto &th : threads) {
    th.join();
  }
  threads.clear();

  std::vector<uint32_t> buffer(n);

  for (int byte_index = 0; byte_index < 4; ++byte_index) {
    const int shift = byte_index * kBitsPerDigit;

    std::vector<std::vector<int>> local_counts(num_threads, std::vector<int>(kRadix, 0));
    for (int t = 0; t < num_threads; ++t) {
      const size_t start = t * chunk_size;
      const size_t end = std::min(start + chunk_size, n);
      threads.emplace_back([&local_counts, &temp, shift, start, end, t]() {
        for (size_t i = start; i < end; ++i) {
          int digit = static_cast<int>((temp[i] >> shift) & 0xFFU);
          ++local_counts[t][digit];
        }
      });
    }
    for (auto &th : threads) {
      th.join();
    }
    threads.clear();

    std::vector<int> global_start(kRadix, 0);
    int current_pos = 0;
    for (int d = 0; d < kRadix; ++d) {
      int total_for_digit = 0;
      for (int t = 0; t < num_threads; ++t) {
        total_for_digit += local_counts[t][d];
      }
      global_start[d] = current_pos;
      current_pos += total_for_digit;
    }

    std::vector<std::vector<int>> thread_offset(num_threads, std::vector<int>(kRadix, 0));
    for (int d = 0; d < kRadix; ++d) {
      int offset = 0;
      for (int t = 0; t < num_threads; ++t) {
        thread_offset[t][d] = offset;
        offset += local_counts[t][d];
      }
    }

    std::vector<std::vector<int>> local_counter(num_threads, std::vector<int>(kRadix, 0));

    for (int t = 0; t < num_threads; ++t) {
      const size_t start = t * chunk_size;
      const size_t end = std::min(start + chunk_size, n);
      threads.emplace_back([&buffer, &temp, &global_start, &thread_offset, &local_counter, shift, start, end, t]() {
        auto &my_counter = local_counter[t];
        for (size_t i = start; i < end; ++i) {
          int digit = static_cast<int>((temp[i] >> shift) & 0xFFU);
          int pos = global_start[digit] + thread_offset[t][digit] + my_counter[digit];
          buffer[static_cast<size_t>(pos)] = temp[i];
          ++my_counter[digit];
        }
      });
    }
    for (auto &th : threads) {
      th.join();
    }
    threads.clear();

    temp.swap(buffer);
  }

  for (int t = 0; t < num_threads; ++t) {
    const size_t start = t * chunk_size;
    const size_t end = std::min(start + chunk_size, n);
    threads.emplace_back([&data, &temp, start, end]() {
      for (size_t i = start; i < end; ++i) {
        data[i] = static_cast<int>(temp[i] ^ kSignMask);
      }
    });
  }
  for (auto &th : threads) {
    th.join();
  }
}

bool ChernovTRadixSortSTL::RunImpl() {
  auto &data = GetOutput();
  if (data.size() <= 1) {
    return true;
  }

  const int num_threads = ppc::util::GetNumThreads();
  RadixSortLSDParallel(data, num_threads);
  return true;
}

bool ChernovTRadixSortSTL::PostProcessingImpl() {
  return std::is_sorted(GetOutput().begin(), GetOutput().end());
}

}  // namespace chernov_t_radix_sort
