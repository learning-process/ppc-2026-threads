#include "sosnina_a_radix_simple_merge/omp/include/ops_omp.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
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

void RadixSortLSD(std::vector<int> &data, std::vector<int> &buffer) {
  size_t idx = 0;
  for (int elem : data) {
    buffer[idx++] = static_cast<int>(static_cast<uint32_t>(elem) ^ kSignFlip);
  }
  std::swap(data, buffer);

  for (int pass = 0; pass < kNumPasses; ++pass) {
    std::vector<int> count(kRadixSize + 1, 0);

    for (auto elem : data) {
      auto digit = static_cast<uint8_t>((static_cast<uint32_t>(elem) >> (pass * kRadixBits)) & 0xFF);
      ++count[digit + 1];
    }

    for (int i = 1; i <= kRadixSize; ++i) {
      count[i] += count[i - 1];
    }

    for (auto elem : data) {
      auto digit = static_cast<uint8_t>((static_cast<uint32_t>(elem) >> (pass * kRadixBits)) & 0xFF);
      buffer[count[digit]++] = elem;
    }

    std::swap(data, buffer);
  }

  for (int &elem : data) {
    elem = static_cast<int>(static_cast<uint32_t>(elem) ^ kSignFlip);
  }
}

void SimpleMerge(const std::vector<int> &left, const std::vector<int> &right, std::vector<int> &result) {
  size_t i = 0;
  size_t j = 0;
  size_t k = 0;

  while (i < left.size() && j < right.size()) {
    if (left[i] <= right[j]) {
      result[k++] = left[i++];
    } else {
      result[k++] = right[j++];
    }
  }

  while (i < left.size()) {
    result[k++] = left[i++];
  }

  while (j < right.size()) {
    result[k++] = right[j++];
  }
}

}  // namespace

SosninaATestTaskOMP::SosninaATestTaskOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = in;
}

bool SosninaATestTaskOMP::ValidationImpl() {
  return !GetInput().empty();
}

bool SosninaATestTaskOMP::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

bool SosninaATestTaskOMP::RunImpl() {
  std::vector<int> &data = GetOutput();
  if (data.size() <= 1) {
    return true;
  }

  const int num_threads = ppc::util::GetNumThreads();
  const int num_parts = std::min(num_threads, static_cast<int>(data.size()));

  if (num_parts <= 1) {
    std::vector<int> buffer(data.size());
    RadixSortLSD(data, buffer);
    return std::ranges::is_sorted(data);
  }

  std::vector<std::vector<int>> parts(num_parts);
  const size_t base_size = data.size() / num_parts;
  const size_t remainder = data.size() % num_parts;
  size_t pos = 0;

  for (int i = 0; i < num_parts; ++i) {
    const size_t part_size = base_size + (std::cmp_less(i, remainder) ? 1 : 0);
    parts[i].assign(data.begin() + static_cast<std::ptrdiff_t>(pos),
                    data.begin() + static_cast<std::ptrdiff_t>(pos + part_size));
    pos += part_size;
  }

#pragma omp parallel for num_threads(num_parts) default(none) shared(parts, num_parts)
  for (int i = 0; i < num_parts; ++i) {
    std::vector<int> buffer(parts[i].size());
    RadixSortLSD(parts[i], buffer);
  }

  std::vector<std::vector<int>> current = std::move(parts);
  while (current.size() > 1) {
    const size_t half = (current.size() + 1) / 2;
    std::vector<std::vector<int>> next(half);

#pragma omp parallel for default(none) shared(current, next) schedule(static)
    for (size_t i = 0; i < current.size() / 2; ++i) {
      next[i].resize(current[2 * i].size() + current[(2 * i) + 1].size());
      SimpleMerge(current[2 * i], current[(2 * i) + 1], next[i]);
    }
    if (current.size() % 2 == 1) {
      next[half - 1] = std::move(current.back());
    }
    current = std::move(next);
  }

  data = std::move(current[0]);
  return std::ranges::is_sorted(data);
}

bool SosninaATestTaskOMP::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace sosnina_a_radix_simple_merge
