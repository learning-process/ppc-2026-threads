#include "sosnina_a_radix_simple_merge/stl/include/ops_stl.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

#include "sosnina_a_radix_simple_merge/common/include/common.hpp"
#include "util/include/util.hpp"

namespace sosnina_a_radix_simple_merge {

namespace {

constexpr int kRadixBits = 8;
constexpr int kRadixSize = 1 << kRadixBits;
constexpr int kNumPasses = sizeof(int) / sizeof(uint8_t);
constexpr uint32_t kSignFlip = 0x80000000U;
constexpr size_t kMinElementsPerPart = 4096;
constexpr size_t kMinElementsPerPartSmall = 32768;
constexpr size_t kSmallArrayThreshold = 1'000'000;
constexpr size_t kLargeArrayThreshold = 20'000'000;

void RadixSortLSD(std::vector<int> &data, std::vector<int> &buffer) {
  size_t idx = 0;
  for (int elem : data) {
    buffer[idx++] = static_cast<int>(static_cast<uint32_t>(elem) ^ kSignFlip);
  }
  std::swap(data, buffer);

  for (int pass = 0; pass < kNumPasses; ++pass) {
    std::array<int, kRadixSize + 1> count{};

    for (auto elem : data) {
      auto digit = static_cast<uint8_t>((static_cast<uint32_t>(elem) >> (pass * kRadixBits)) & 0xFF);
      ++count.at(static_cast<size_t>(digit) + 1U);
    }

    for (int i = 1; i <= kRadixSize; ++i) {
      const auto ui = static_cast<size_t>(i);
      count.at(ui) += count.at(ui - 1U);
    }

    for (auto elem : data) {
      auto digit = static_cast<uint8_t>((static_cast<uint32_t>(elem) >> (pass * kRadixBits)) & 0xFF);
      const auto di = static_cast<size_t>(digit);
      const int write_pos = count.at(di)++;
      buffer[static_cast<size_t>(write_pos)] = elem;
    }

    std::swap(data, buffer);
  }

  for (int &elem : data) {
    elem = static_cast<int>(static_cast<uint32_t>(elem) ^ kSignFlip);
  }
}

void SimpleMerge(const std::vector<int> &left, const std::vector<int> &right, std::vector<int> &result) {
  std::ranges::merge(left, right, result.begin());
}

/// Диапазон [begin, end) по индексам, до num_threads потоков (как грубый аналог parallel for).
template <typename F>
void ParallelForRange(size_t begin, size_t end, int num_threads, F &&fn) {
  if (begin >= end) {
    return;
  }
  std::decay_t<F> func{std::forward<F>(fn)};
  num_threads = std::max(1, std::min(num_threads, static_cast<int>(end - begin)));
  const size_t n = end - begin;
  const size_t chunk = (n + static_cast<size_t>(num_threads) - 1) / static_cast<size_t>(num_threads);
  std::vector<std::thread> threads;
  for (int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
    const size_t lo = begin + (static_cast<size_t>(thread_idx) * chunk);
    if (lo >= end) {
      break;
    }
    const size_t hi = std::min(end, lo + chunk);
    threads.emplace_back([lo, hi, &func]() {
      for (size_t i = lo; i < hi; ++i) {
        func(i);
      }
    });
  }
  for (auto &th : threads) {
    th.join();
  }
}

}  // namespace

SosninaATestTaskSTL::SosninaATestTaskSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = in;
}

bool SosninaATestTaskSTL::ValidationImpl() {
  return !GetInput().empty();
}

bool SosninaATestTaskSTL::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool SosninaATestTaskSTL::RunImpl() {
  std::vector<int> &data = GetOutput();
  if (data.size() <= 1) {
    return true;
  }

  const int num_threads = ppc::util::GetNumThreads();
  const size_t per_thread_floor = data.size() >= kLargeArrayThreshold
                                      ? (data.size() / static_cast<size_t>(std::max(1, num_threads)))
                                      : (data.size() / static_cast<size_t>(std::max(1, 2 * num_threads)));
  size_t min_chunk_base = kMinElementsPerPart;
  if (data.size() < kMinElementsPerPartSmall) {
    min_chunk_base = std::max(size_t{1}, data.size() / static_cast<size_t>(std::max(1, 2 * num_threads)));
  } else if (data.size() < kSmallArrayThreshold) {
    min_chunk_base = kMinElementsPerPartSmall;
  }
  const size_t min_chunk = std::max(min_chunk_base, per_thread_floor);
  const int max_parts_by_grain = std::max(1, static_cast<int>(data.size() / min_chunk));
  const int num_parts = std::min({num_threads, static_cast<int>(data.size()), max_parts_by_grain});

  if (num_parts <= 1) {
    std::vector<int> buffer(data.size());
    RadixSortLSD(data, buffer);
    return std::ranges::is_sorted(data);
  }

  std::vector<std::vector<int>> parts(static_cast<size_t>(num_parts));
  const size_t base_size = data.size() / static_cast<size_t>(num_parts);
  const size_t remainder = data.size() % static_cast<size_t>(num_parts);
  size_t pos = 0;

  for (int i = 0; i < num_parts; ++i) {
    const size_t part_size = base_size + (std::cmp_less(i, remainder) ? 1U : 0U);
    parts[static_cast<size_t>(i)].assign(data.begin() + static_cast<std::ptrdiff_t>(pos),
                                         data.begin() + static_cast<std::ptrdiff_t>(pos + part_size));
    pos += part_size;
  }

  ParallelForRange(0, static_cast<size_t>(num_parts), num_threads, [&](size_t i) {
    auto &part = parts[i];
    std::vector<int> buffer(part.size());
    RadixSortLSD(part, buffer);
  });

  std::vector<std::vector<int>> current = std::move(parts);
  while (current.size() > 1) {
    const size_t half = (current.size() + 1) / 2;
    std::vector<std::vector<int>> next(half);

    const size_t pair_count = current.size() / 2;
    ParallelForRange(0, pair_count, num_threads, [&](size_t idx) {
      auto &left = current[2 * idx];
      auto &right = current[(2 * idx) + 1];
      next[idx].resize(left.size() + right.size());
      SimpleMerge(left, right, next[idx]);
      std::vector<int>().swap(left);
      std::vector<int>().swap(right);
    });
    if (current.size() % 2 == 1) {
      next[half - 1] = std::move(current.back());
    }
    current = std::move(next);
  }

  data = std::move(current[0]);
  return std::ranges::is_sorted(data);
}

bool SosninaATestTaskSTL::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace sosnina_a_radix_simple_merge
