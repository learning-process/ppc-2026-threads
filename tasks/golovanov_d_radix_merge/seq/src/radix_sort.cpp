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

    if ((bits & kSignMask) != 0U) {  // отрицательное
      bits = ~bits;
    } else {  // положительное
      bits ^= kSignMask;
    }

    data[i] = bits;
  }

  std::vector<uint64_t> buffer(n);

  for (int byte = 0; byte < kBytes; ++byte) {
    std::array<size_t, kRadix> count{};

    for (size_t i = 0; i < n; ++i) {
      size_t b = static_cast<size_t>((data[i] >> (byte * 8)) & kByteMask);
      count.at(b)++;
    }

    size_t sum = 0;
    for (auto &c : count) {
      size_t tmp = c;
      c = sum;
      sum += tmp;
    }

    for (size_t i = 0; i < n; ++i) {
      size_t b = static_cast<size_t>((data[i] >> (byte * 8)) & kByteMask);
      buffer[count.at(b)++] = data[i];
    }

    data.swap(buffer);
  }

  for (size_t i = 0; i < n; ++i) {
    uint64_t bits = data[i];

    if ((bits & kSignMask) != 0U) {
      bits ^= kSignMask;
    } else {
      bits = ~bits;
    }

    std::memcpy(&arr[i], &bits, sizeof(double));
  }
}
