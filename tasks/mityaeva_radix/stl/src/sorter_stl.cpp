#include "mityaeva_radix/stl/include/sorter_stl.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <thread>
#include <vector>

#include "util/include/util.hpp"

namespace mityaeva_radix {

namespace {
template <typename Index, typename Functor>
void ParallelForByThreads(Index start, Index end, Functor f, size_t num_threads) {
  if (num_threads < 2) {
    for (auto i = start; i < end; i++) {
      f(i);
    }
    return;
  }

  std::vector<std::thread> threads;
  threads.reserve(num_threads);
  for (auto thread_idx = start; thread_idx < end; thread_idx++) {
    threads.emplace_back([&, thread_idx] { f(thread_idx); });
  }

  for (auto &thread : threads) {
    thread.join();
  }
}

template <typename Index, typename Functor>
void ParallelForByData(Index start, Index end, Functor f, size_t num_threads) {
  if (num_threads < 2) {
    f(start, end);
    return;
  }
  auto portion = (end - start) / num_threads;
  auto remainder = (end - start) % num_threads;

  std::vector<std::thread> threads;
  threads.reserve(num_threads);
  for (auto thread_idx = 0UZ; thread_idx < num_threads; thread_idx++) {
    auto begin = start + (portion * thread_idx) + std::min(thread_idx, remainder);
    auto end = start + ((thread_idx + 1) * portion) + std::min(thread_idx + 1, remainder);
    threads.emplace_back([&, begin, end] { f(begin, end); });
  }

  for (auto &thread : threads) {
    thread.join();
  }
}
}  // namespace

uint64_t SorterStl::DoubleToSortable(uint64_t x) {
  if ((x & 0x8000000000000000ULL) != 0U) {
    return ~x;
  }
  return x | 0x8000000000000000ULL;
}

uint64_t SorterStl::SortableToDouble(uint64_t x) {
  if ((x & 0x8000000000000000ULL) != 0U) {
    return x & 0x7FFFFFFFFFFFFFFFULL;
  }
  return ~x;
}

void SorterStl::CountingPass(std::vector<uint64_t> *current, std::vector<uint64_t> *next, int shift, int radix,
                             int num_threads, size_t data_size) {
  std::vector<std::vector<int>> thread_counters(num_threads, std::vector<int>(radix, 0));
  ParallelForByThreads(0, num_threads, [&](int thread_id) {
    size_t chunk_size = data_size / static_cast<size_t>(num_threads);
    size_t start = static_cast<size_t>(thread_id) * chunk_size;
    size_t end = (thread_id == num_threads - 1) ? data_size : start + chunk_size;

    auto &local_counters = thread_counters[thread_id];

    for (size_t i = start; i < end; i++) {
      int digit = static_cast<int>(((*current)[i] >> static_cast<size_t>(shift)) & static_cast<size_t>(radix - 1));
      local_counters[digit]++;
    }
  }, num_threads);

  std::vector<int> prefix_sums(static_cast<size_t>(radix * num_threads), 0);

  int total = 0;
  for (int digit = 0; digit < radix; digit++) {
    int digit_sum = 0;
    for (int thread_idx = 0; thread_idx < num_threads; thread_idx++) {
      prefix_sums[(thread_idx * radix) + digit] = total + digit_sum;
      digit_sum += thread_counters[thread_idx][digit];
    }
    total += digit_sum;
  }

  ParallelForByThreads(0, num_threads, [&](int thread_id) {
    size_t chunk_size = data_size / static_cast<size_t>(num_threads);
    size_t start = static_cast<size_t>(thread_id) * chunk_size;
    size_t end = (thread_id == num_threads - 1) ? data_size : start + chunk_size;

    std::vector<int> local_pos(radix, 0);
    for (int digit = 0; digit < radix; digit++) {
      local_pos[digit] = prefix_sums[(thread_id * radix) + digit];
    }

    for (size_t i = start; i < end; i++) {
      int digit = static_cast<int>(((*current)[i] >> static_cast<size_t>(shift)) & static_cast<size_t>(radix - 1));
      auto pos = static_cast<size_t>(local_pos[digit]++);
      (*next)[pos] = (*current)[i];
    }
  }, num_threads);
}

void SorterStl::Sort(std::vector<double> &data) {
  if (data.size() <= 1) {
    return;
  }
  int num_threads = ppc::util::GetNumThreads();
  std::vector<double> temp(data.size());
  std::vector<uint64_t> as_uint(data.size());

  ParallelForByData(0UZ, data.size(), [&](size_t start, size_t end) {
    for (auto i = start; i < end; i++) {
      uint64_t bits = 0;
      std::memcpy(&bits, &data[i], sizeof(double));
      as_uint[i] = DoubleToSortable(bits);
    }
  }, num_threads);
  const int bits_per_pass = 8;
  const int radix = 1 << bits_per_pass;
  const int passes = static_cast<int>(sizeof(uint64_t) * 8 / bits_per_pass);
  std::vector<uint64_t> uint_temp(data.size());
  std::vector<uint64_t> *current = &as_uint;
  std::vector<uint64_t> *next = &uint_temp;

  for (int pass = 0; pass < passes; pass++) {
    int shift = pass * bits_per_pass;
    CountingPass(current, next, shift, radix, num_threads, data.size());
    std::swap(current, next);
  }

  ParallelForByData(0UZ, data.size(), [&](size_t start, size_t end) {
    for (auto i = start; i < end; i++) {
      uint64_t bits = 0;
      if (current == &as_uint) {
        bits = SortableToDouble(as_uint[i]);
      } else {
        bits = SortableToDouble(uint_temp[i]);
      }
      std::memcpy(&data[i], &bits, sizeof(double));
    }
  }, num_threads);
}

}  // namespace mityaeva_radix
