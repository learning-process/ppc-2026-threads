#include "mityaeva_radix/tbb/include/sorter_tbb.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>

#include "oneapi/tbb/parallel_for.h"
#include "oneapi/tbb/partitioner.h"
#include "util/include/util.hpp"

namespace mityaeva_radix {

uint64_t SorterTbb::DoubleToSortable(uint64_t x) {
  if ((x & 0x8000000000000000ULL) != 0U) {
    return ~x;
  }
  return x | 0x8000000000000000ULL;
}

uint64_t SorterTbb::SortableToDouble(uint64_t x) {
  if ((x & 0x8000000000000000ULL) != 0U) {
    return x & 0x7FFFFFFFFFFFFFFFULL;
  }
  return ~x;
}

void SorterTbb::CountingPass(std::vector<uint64_t> *current, std::vector<uint64_t> *next, int shift, int radix,
                             int num_threads, size_t data_size) {
  std::vector<std::vector<int>> thread_counters(num_threads, std::vector<int>(radix, 0));

  tbb::parallel_for(0, num_threads, [&](int thread_id) {
    size_t chunk_size = data_size / static_cast<size_t>(num_threads);
    size_t start = static_cast<size_t>(thread_id) * chunk_size;
    size_t end = (thread_id == num_threads - 1) ? data_size : start + chunk_size;

    auto &local_counters = thread_counters[thread_id];

    for (size_t i = start; i < end; i++) {
      int digit = static_cast<int>(((*current)[i] >> static_cast<size_t>(shift)) & static_cast<size_t>(radix - 1));
      local_counters[digit]++;
    }
  }, tbb::static_partitioner());

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

  tbb::parallel_for(0, num_threads, [&](int thread_id) {
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
  }, tbb::static_partitioner());
}

void SorterTbb::Sort(std::vector<double> &data) {
  if (data.size() <= 1) {
    return;
  }

  int num_threads = ppc::util::GetNumThreads();
  std::vector<double> temp(data.size());
  std::vector<uint64_t> as_uint(data.size());

  tbb::parallel_for(tbb::blocked_range<std::size_t>(0, data.size()),
                    [&data, &as_uint](const tbb::blocked_range<std::size_t> &range) {
    for (auto i = range.begin(); i < range.end(); i++) {
      uint64_t bits = 0;
      std::memcpy(&bits, &data[i], sizeof(double));
      as_uint[i] = DoubleToSortable(bits);
    }
  });
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

  tbb::parallel_for(tbb::blocked_range<size_t>(0, data.size()),
                    [&data, &current, &as_uint, &uint_temp](const tbb::blocked_range<std::size_t> &range) {
    for (auto i = range.begin(); i < range.end(); i++) {
      uint64_t bits = 0;
      if (current == &as_uint) {
        bits = SortableToDouble(as_uint[i]);
      } else {
        bits = SortableToDouble(uint_temp[i]);
      }
      std::memcpy(&data[i], &bits, sizeof(double));
    }
  });
}

}  // namespace mityaeva_radix
