#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <utility>
#include <vector>

namespace egashin_k_radix_simple_merge::radix_utils {

namespace detail {

constexpr uint64_t kSignBit = 0x8000000000000000ULL;

inline uint64_t ToSortable(double value) {
  uint64_t bits = 0;
  std::memcpy(&bits, &value, sizeof(value));
  return ((bits & kSignBit) != 0U) ? ~bits : (bits ^ kSignBit);
}

inline double FromSortable(uint64_t key) {
  const uint64_t bits = ((key & kSignBit) != 0U) ? (key ^ kSignBit) : ~key;
  double value = 0.0;
  std::memcpy(&value, &bits, sizeof(value));
  return value;
}

inline void CountingPass(const std::vector<uint64_t> &source, std::vector<uint64_t> &destination, int byte_index) {
  std::array<size_t, 256> count{};

  for (uint64_t value : source) {
    const auto byte = static_cast<uint8_t>((value >> (byte_index * 8)) & 0xFFU);
    count.at(byte)++;
  }

  std::array<size_t, 256> position{};
  for (size_t i = 1; i < count.size(); ++i) {
    position.at(i) = position.at(i - 1) + count.at(i - 1);
  }

  for (uint64_t value : source) {
    const auto byte = static_cast<uint8_t>((value >> (byte_index * 8)) & 0xFFU);
    destination[position.at(byte)++] = value;
  }
}

}  // namespace detail

inline int WorkerCount(size_t size, int requested) {
  if (size == 0) {
    return 1;
  }
  return std::min(std::max(requested, 1), static_cast<int>(size));
}

inline std::vector<std::pair<size_t, size_t>> MakeRanges(size_t size, int workers) {
  std::vector<std::pair<size_t, size_t>> ranges(static_cast<size_t>(workers));
  const auto worker_count = static_cast<size_t>(workers);
  const size_t base = size / worker_count;
  const size_t extra = size % worker_count;

  size_t begin = 0;
  for (size_t i = 0; i < worker_count; ++i) {
    const size_t block = base + (i < extra ? 1 : 0);
    ranges[i] = {begin, begin + block};
    begin += block;
  }

  return ranges;
}

inline void SortRange(std::vector<double> &data, size_t left, size_t right) {
  if (right - left < 2) {
    return;
  }

  const size_t size = right - left;
  std::vector<uint64_t> keys(size);
  std::vector<uint64_t> buffer(size);

  for (size_t i = 0; i < size; ++i) {
    keys[i] = detail::ToSortable(data[left + i]);
  }

  auto *source = &keys;
  auto *destination = &buffer;
  for (int byte_index = 0; byte_index < 8; ++byte_index) {
    detail::CountingPass(*source, *destination, byte_index);
    std::swap(source, destination);
  }

  for (size_t i = 0; i < size; ++i) {
    data[left + i] = detail::FromSortable((*source)[i]);
  }
}

inline std::vector<double> Merge(const std::vector<double> &left, const std::vector<double> &right) {
  std::vector<double> result;
  result.reserve(left.size() + right.size());

  size_t i = 0;
  size_t j = 0;
  while (i < left.size() && j < right.size()) {
    if (left[i] <= right[j]) {
      result.push_back(left[i++]);
    } else {
      result.push_back(right[j++]);
    }
  }

  result.insert(result.end(), left.begin() + static_cast<std::ptrdiff_t>(i), left.end());
  result.insert(result.end(), right.begin() + static_cast<std::ptrdiff_t>(j), right.end());
  return result;
}

inline std::vector<std::vector<double>> MakeParts(const std::vector<double> &data,
                                                  const std::vector<std::pair<size_t, size_t>> &ranges) {
  std::vector<std::vector<double>> parts(ranges.size());
  for (size_t i = 0; i < ranges.size(); ++i) {
    parts[i] = std::vector<double>(data.begin() + static_cast<std::ptrdiff_t>(ranges[i].first),
                                   data.begin() + static_cast<std::ptrdiff_t>(ranges[i].second));
  }
  return parts;
}

}  // namespace egashin_k_radix_simple_merge::radix_utils
