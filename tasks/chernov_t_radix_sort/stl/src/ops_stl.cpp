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

void ChernovTRadixSortSTL::RadixSortLSDSequential(std::vector<int> &data) {
  if (data.size() <= 1) {
    return;
  }

  const size_t n = data.size();
  std::vector<uint32_t> temp(n);
  for (size_t i = 0; i < n; ++i) {
    temp[i] = static_cast<uint32_t>(data[i]) ^ kSignMask;
  }

  std::vector<uint32_t> buffer(n);

  for (int byte = 0; byte < 4; ++byte) {
    const int shift = byte * kBitsPerDigit;
    std::vector<int> count(kRadix, 0);

    for (size_t i = 0; i < n; ++i) {
      ++count[static_cast<int>((temp[i] >> shift) & 0xFFU)];
    }
    for (int i = 1; i < kRadix; ++i) {
      count[i] += count[i - 1];
    }
    for (size_t i = n; i-- > 0;) {
      buffer[static_cast<size_t>(--count[static_cast<int>((temp[i] >> shift) & 0xFFU)])] = temp[i];
    }
    temp.swap(buffer);
  }

  for (size_t i = 0; i < n; ++i) {
    data[i] = static_cast<int>(temp[i] ^ kSignMask);
  }
}

void ChernovTRadixSortSTL::RadixSortLSDParallel(std::vector<int> &data, int num_threads) {
  const size_t n = data.size();
  if (n <= 1) {
    return;
  }

  std::vector<uint32_t> temp(n);
  std::vector<std::thread> threads;
  const size_t chunk_size = (n + static_cast<size_t>(num_threads) - 1) / static_cast<size_t>(num_threads);

  for (int t = 0; t < num_threads; ++t) {
    const size_t start = static_cast<size_t>(t) * chunk_size;
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

    std::vector<std::vector<int>> local_counts(static_cast<size_t>(num_threads), std::vector<int>(kRadix, 0));
    for (int t = 0; t < num_threads; ++t) {
      const size_t start = static_cast<size_t>(t) * chunk_size;
      const size_t end = std::min(start + chunk_size, n);
      threads.emplace_back([t, &local_counts, &temp, shift, start, end]() {
        auto &cnt = local_counts[static_cast<size_t>(t)];
        for (size_t i = start; i < end; ++i) {
          ++cnt[static_cast<int>((temp[i] >> shift) & 0xFFU)];
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
        total_for_digit += local_counts[static_cast<size_t>(t)][d];
      }
      global_start[d] = current_pos;
      current_pos += total_for_digit;
    }

    std::vector<std::vector<int>> thread_offset(static_cast<size_t>(num_threads), std::vector<int>(kRadix, 0));
    for (int d = 0; d < kRadix; ++d) {
      int offset = 0;
      for (int t = 0; t < num_threads; ++t) {
        thread_offset[static_cast<size_t>(t)][d] = offset;
        offset += local_counts[static_cast<size_t>(t)][d];
      }
    }

    std::vector<std::vector<int>> local_counter(static_cast<size_t>(num_threads), std::vector<int>(kRadix, 0));

    for (int t = 0; t < num_threads; ++t) {
      const size_t start = static_cast<size_t>(t) * chunk_size;
      const size_t end = std::min(start + chunk_size, n);
      threads.emplace_back([t, &buffer, &temp, &global_start, &thread_offset, &local_counter, shift, start, end]() {
        auto &my_counter = local_counter[static_cast<size_t>(t)];
        for (size_t i = start; i < end; ++i) {
          int digit = static_cast<int>((temp[i] >> shift) & 0xFFU);
          int pos = global_start[digit] + thread_offset[static_cast<size_t>(t)][digit] + my_counter[digit];
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
    const size_t start = static_cast<size_t>(t) * chunk_size;
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

  if (data.size() < 1000) {
    RadixSortLSDSequential(data);
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
