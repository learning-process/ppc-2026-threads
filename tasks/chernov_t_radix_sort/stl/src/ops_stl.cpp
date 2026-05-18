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
constexpr int kRadix = 1 << kBitsPerDigit;
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

void ChernovTRadixSortSTL::ConvertToIntegers(const std::vector<int> &input, std::vector<uint32_t> &output, size_t start,
                                             size_t end, uint32_t sign_mask) {
  for (size_t i = start; i < end; ++i) {
    output[i] = static_cast<uint32_t>(input[i]) ^ sign_mask;
  }
}

void ChernovTRadixSortSTL::ComputeLocalHistograms(const std::vector<uint32_t> &data,
                                                  std::vector<std::vector<int>> &local_counts, size_t start, size_t end,
                                                  int shift, int thread_idx) {
  auto &cnt = local_counts[static_cast<size_t>(thread_idx)];
  for (size_t i = start; i < end; ++i) {
    int digit = static_cast<int>((data[i] >> shift) & 0xFFU);
    ++cnt[digit];
  }
}

void ChernovTRadixSortSTL::ComputeGlobalStarts(const std::vector<std::vector<int>> &local_counts,
                                               std::vector<int> &global_start, int k_radix, int num_threads) {
  int current_pos = 0;
  for (int digit_idx = 0; digit_idx < k_radix; ++digit_idx) {
    int total = 0;
    for (int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
      total += local_counts[static_cast<size_t>(thread_idx)][digit_idx];
    }
    global_start[digit_idx] = current_pos;
    current_pos += total;
  }
}

void ChernovTRadixSortSTL::ComputeThreadOffsets(const std::vector<std::vector<int>> &local_counts,
                                                std::vector<std::vector<int>> &thread_offset, int k_radix,
                                                int num_threads) {
  for (int digit_idx = 0; digit_idx < k_radix; ++digit_idx) {
    int offset = 0;
    for (int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
      thread_offset[static_cast<size_t>(thread_idx)][digit_idx] = offset;
      offset += local_counts[static_cast<size_t>(thread_idx)][digit_idx];
    }
  }
}

void ChernovTRadixSortSTL::ScatterElements(const std::vector<uint32_t> &input, std::vector<uint32_t> &output,
                                           const std::vector<int> &global_start,
                                           const std::vector<std::vector<int>> &thread_offset,
                                           std::vector<std::vector<int>> &local_counter, size_t start, size_t end,
                                           int shift, int thread_idx) {
  auto &my_counter = local_counter[static_cast<size_t>(thread_idx)];
  for (size_t i = start; i < end; ++i) {
    int digit = static_cast<int>((input[i] >> shift) & 0xFFU);
    int pos = global_start[digit] + thread_offset[static_cast<size_t>(thread_idx)][digit] + my_counter[digit];
    output[static_cast<size_t>(pos)] = input[i];
    ++my_counter[digit];
  }
}

void ChernovTRadixSortSTL::ConvertFromIntegers(const std::vector<uint32_t> &input, std::vector<int> &output,
                                               size_t start, size_t end, uint32_t sign_mask) {
  for (size_t i = start; i < end; ++i) {
    output[i] = static_cast<int>(input[i] ^ sign_mask);
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

  for (int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    const size_t start = static_cast<size_t>(thread_idx) * chunk_size;
    const size_t end = std::min(start + chunk_size, n);
    threads.emplace_back([&data, &temp, start, end]() { ConvertToIntegers(data, temp, start, end, kSignMask); });
  }
  for (auto &th : threads) {
    th.join();
  }
  threads.clear();

  std::vector<uint32_t> buffer(n);

  for (int byte_index = 0; byte_index < 4; ++byte_index) {
    const int shift = byte_index * kBitsPerDigit;

    // Гистограммы
    std::vector<std::vector<int>> local_counts(static_cast<size_t>(num_threads), std::vector<int>(kRadix, 0));
    for (int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
      const size_t start = static_cast<size_t>(thread_idx) * chunk_size;
      const size_t end = std::min(start + chunk_size, n);
      threads.emplace_back([&temp, &local_counts, start, end, shift, thread_idx]() {
        ComputeLocalHistograms(temp, local_counts, start, end, shift, thread_idx);
      });
    }
    for (auto &th : threads) {
      th.join();
    }
    threads.clear();

    std::vector<int> global_start(kRadix, 0);
    ComputeGlobalStarts(local_counts, global_start, kRadix, num_threads);

    std::vector<std::vector<int>> thread_offset(static_cast<size_t>(num_threads), std::vector<int>(kRadix, 0));
    ComputeThreadOffsets(local_counts, thread_offset, kRadix, num_threads);

    std::vector<std::vector<int>> local_counter(static_cast<size_t>(num_threads), std::vector<int>(kRadix, 0));

    for (int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
      const size_t start = static_cast<size_t>(thread_idx) * chunk_size;
      const size_t end = std::min(start + chunk_size, n);
      threads.emplace_back(
          [&temp, &buffer, &global_start, &thread_offset, &local_counter, start, end, shift, thread_idx]() {
        ScatterElements(temp, buffer, global_start, thread_offset, local_counter, start, end, shift, thread_idx);
      });
    }
    for (auto &th : threads) {
      th.join();
    }
    threads.clear();

    temp.swap(buffer);
  }

  for (int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    const size_t start = static_cast<size_t>(thread_idx) * chunk_size;
    const size_t end = std::min(start + chunk_size, n);
    threads.emplace_back([&temp, &data, start, end]() { ConvertFromIntegers(temp, data, start, end, kSignMask); });
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
