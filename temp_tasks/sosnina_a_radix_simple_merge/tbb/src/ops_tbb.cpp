#include "sosnina_a_radix_simple_merge/tbb/include/ops_tbb.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "oneapi/tbb/parallel_for.h"
#include "oneapi/tbb/partitioner.h"
#include "sosnina_a_radix_simple_merge/common/include/common.hpp"
#include "util/include/util.hpp"

namespace sosnina_a_radix_simple_merge {

namespace {

constexpr int kRadixBits = 8;
constexpr int kRadixSize = 1 << kRadixBits;
constexpr int kNumPasses = sizeof(int) / sizeof(uint8_t);
constexpr uint32_t kSignFlip = 0x80000000U;
constexpr size_t kMinElementsPerPart = 4096;
/// На n < порога — крупнее части, меньше уровней merge (лучше E на малых входах).
constexpr size_t kMinElementsPerPartSmall = 32768;
constexpr size_t kSmallArrayThreshold = 1'000'000;
/// От этого размера — «крупный» вход: целимся в ~не больше одной части на поток, глубже merge дороже.
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

}  // namespace

SosninaATestTaskTBB::SosninaATestTaskTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = in;
}

bool SosninaATestTaskTBB::ValidationImpl() {
  return !GetInput().empty();
}

bool SosninaATestTaskTBB::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool SosninaATestTaskTBB::RunImpl() {
  std::vector<int> &data = GetOutput();
  if (data.size() <= 1) {
    return true;
  }

  const int num_threads = ppc::util::GetNumThreads();
  const size_t min_chunk_base = (data.size() < kSmallArrayThreshold) ? kMinElementsPerPartSmall : kMinElementsPerPart;
  // На больших массивах — не мельче ~ n/T (меньше лишних уровней merge, толще radix-куски).
  // На малых — оставляем запас n/(2T), чтобы не раздувать число частей.
  const size_t per_thread_floor = data.size() >= kLargeArrayThreshold
                                      ? (data.size() / static_cast<size_t>(std::max(1, num_threads)))
                                      : (data.size() / static_cast<size_t>(std::max(1, 2 * num_threads)));
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

  tbb::parallel_for(0, num_parts, [&](int i) {
    std::vector<int> buffer(parts[static_cast<size_t>(i)].size());
    RadixSortLSD(parts[static_cast<size_t>(i)], buffer);
  }, tbb::simple_partitioner{});

  std::vector<std::vector<int>> current = std::move(parts);
  while (current.size() > 1) {
    const size_t half = (current.size() + 1) / 2;
    std::vector<std::vector<int>> next(half);

    const size_t pair_count = current.size() / 2;
    tbb::parallel_for(static_cast<size_t>(0), pair_count, [&](size_t idx) {
      std::vector<int> &left = current[2 * idx];
      std::vector<int> &right = current[(2 * idx) + 1];
      next[idx].resize(left.size() + right.size());
      SimpleMerge(left, right, next[idx]);
      std::vector<int>().swap(left);
      std::vector<int>().swap(right);
    }, tbb::simple_partitioner{});
    if (current.size() % 2 == 1) {
      next[half - 1] = std::move(current.back());
    }
    current = std::move(next);
  }

  data = std::move(current[0]);
  return std::ranges::is_sorted(data);
}

bool SosninaATestTaskTBB::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace sosnina_a_radix_simple_merge
