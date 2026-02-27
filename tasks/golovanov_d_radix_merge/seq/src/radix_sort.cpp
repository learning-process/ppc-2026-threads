#include "../include/radix_sort.hpp"

#include <array>
#include <cstdint>
#include <cstring>
#include <vector>

namespace {
constexpr int kBytes = 8;
constexpr std::size_t kRadix = 256;
constexpr std::uint64_t kSignMask = 1ULL << 63;
constexpr std::uint64_t kByteMask = 0xFFULL;
uint64_t ToSortable(uint64_t bits, const uint64_t sign_mask) {
  return ((bits & sign_mask) != 0U) ? ~bits : (bits ^ sign_mask);
}

uint64_t FromSortable(uint64_t bits, const uint64_t sign_mask) {
  return ((bits & sign_mask) != 0U) ? (bits ^ sign_mask) : ~bits;
}
}  // namespace

void RadixSort::Sort(std::vector<double> &arr) {
  const size_t n = arr.size();
  if (n == 0) {
    return;
  }

  std::vector<uint64_t> data(n);

  for (size_t i = 0; i < n; ++i) {
    uint64_t bits = 0;
    std::memcpy(&bits, &arr[i], sizeof(double));

    bits = ToSortable(bits, kSignMask);

    data[i] = bits;
  }

  std::vector<uint64_t> buffer(n);

  for (int byte = 0; byte < kBytes; ++byte) {
    std::array<size_t, kRadix> count{};

    for (size_t i = 0; i < n; ++i) {
      const auto b = static_cast<size_t>((data[i] >> (byte * 8)) & kByteMask);
      count.at(b)++;
    }

    size_t sum = 0;
    for (auto &c : count) {
      size_t tmp = c;
      c = sum;
      sum += tmp;
    }

    for (size_t i = 0; i < n; ++i) {
      const auto b = static_cast<size_t>((data[i] >> (byte * 8)) & kByteMask);
      buffer[count.at(b)++] = data[i];
    }

    data.swap(buffer);
  }

  for (size_t i = 0; i < n; ++i) {
    uint64_t bits = data[i];

    bits = FromSortable(bits, kSignMask);

    std::memcpy(&arr[i], &bits, sizeof(double));
  }
}
