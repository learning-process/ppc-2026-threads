#include "frolova_s_radix_sort_double/tbb/include/ops_tbb.hpp"

#include <oneapi/tbb/parallel_for.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "frolova_s_radix_sort_double/common/include/common.hpp"

namespace frolova_s_radix_sort_double {

FrolovaSRadixSortDoubleTBB::FrolovaSRadixSortDoubleTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool FrolovaSRadixSortDoubleTBB::ValidationImpl() {
  return !GetInput().empty();
}

bool FrolovaSRadixSortDoubleTBB::PreProcessingImpl() {
  return true;
}

bool FrolovaSRadixSortDoubleTBB::RunImpl() {
  const std::vector<double> &input = GetInput();
  if (input.empty()) {
    return false;
  }

  std::vector<double> working = input;
  const std::size_t n = working.size();

  constexpr int kRadix = 256;
  constexpr int kNumBits = 8;
  constexpr int kNumPasses = sizeof(std::uint64_t);

  std::vector<double> temp(n);

  for (int pass = 0; pass < kNumPasses; ++pass) {
    std::array<std::atomic<int>, kRadix> count{};
    tbb::parallel_for(std::size_t{0}, n, [&](std::size_t i) {
      auto bits = std::bit_cast<std::uint64_t>(working[i]);
      int byte = static_cast<int>((bits >> (pass * kNumBits)) & 0xFF);
      count.at(byte).fetch_add(1, std::memory_order_relaxed);
    });

    std::array<int, kRadix> offset{};
    int total = 0;
    for (int i = 0; i < kRadix; ++i) {
      offset.at(i) = total;
      total += count.at(i).load();
    }

    for (double val : working) {
      auto bits = std::bit_cast<std::uint64_t>(val);
      int byte = static_cast<int>((bits >> (pass * kNumBits)) & 0xFF);
      temp.at(offset.at(byte)++) = val;
    }

    working.swap(temp);
  }

  std::vector<double> negative;
  std::vector<double> positive;
  negative.reserve(n);
  positive.reserve(n);

  for (double val : working) {
    if (val < 0.0) {
      negative.push_back(val);
    } else {
      positive.push_back(val);
    }
  }
  std::ranges::reverse(negative);

  working.clear();
  working.insert(working.end(), negative.begin(), negative.end());
  working.insert(working.end(), positive.begin(), positive.end());

  GetOutput() = std::move(working);
  return true;
}

bool FrolovaSRadixSortDoubleTBB::PostProcessingImpl() {
  return true;
}

}  // namespace frolova_s_radix_sort_double
