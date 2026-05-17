#include "frolova_s_radix_sort_double/stl/include/ops_stl.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <thread>
#include <utility>
#include <vector>

#include "frolova_s_radix_sort_double/common/include/common.hpp"

namespace frolova_s_radix_sort_double {

namespace {

constexpr int kRadix = 256;
constexpr int kNumBits = 8;
constexpr int kNumPasses = sizeof(std::uint64_t);

std::array<int, kRadix> ComputeHistogram(const std::vector<double> &data, int pass, unsigned int num_threads) {
  std::array<std::atomic<int>, kRadix> count{};
  const std::size_t n = data.size();
  std::vector<std::thread> threads;
  std::size_t chunk_size = (n + num_threads - 1) / num_threads;

  for (unsigned int tid = 0; tid < num_threads; ++tid) {
    std::size_t start = tid * chunk_size;
    if (start >= n) {
      break;
    }
    std::size_t end = std::min(start + chunk_size, n);
    threads.emplace_back([&, start, end]() {
      for (std::size_t i = start; i < end; ++i) {
        auto bits = std::bit_cast<std::uint64_t>(data[i]);
        int byte = static_cast<int>((bits >> (pass * kNumBits)) & 0xFF);
        count.at(byte).fetch_add(1, std::memory_order_relaxed);
      }
    });
  }
  for (auto &t : threads) {
    t.join();
  }

  std::array<int, kRadix> result{};
  for (int i = 0; i < kRadix; ++i) {
    result.at(i) = count.at(i).load();
  }
  return result;
}

std::array<int, kRadix> BuildOffsets(const std::array<int, kRadix> &histogram) {
  std::array<int, kRadix> offsets{};
  int total = 0;
  for (int i = 0; i < kRadix; ++i) {
    offsets.at(i) = total;
    total += histogram.at(i);
  }
  return offsets;
}

void Distribute(const std::vector<double> &src, std::vector<double> &dst, std::array<int, kRadix> &offsets, int pass) {
  for (double val : src) {
    auto bits = std::bit_cast<std::uint64_t>(val);
    int byte = static_cast<int>((bits >> (pass * kNumBits)) & 0xFF);
    dst.at(offsets.at(byte)++) = val;
  }
}

void FixNegativeOrder(std::vector<double> &data) {
  std::vector<double> negative;
  std::vector<double> positive;
  negative.reserve(data.size());
  positive.reserve(data.size());

  for (double val : data) {
    if (val < 0.0) {
      negative.push_back(val);
    } else {
      positive.push_back(val);
    }
  }
  std::ranges::reverse(negative);

  data.clear();
  data.insert(data.end(), negative.begin(), negative.end());
  data.insert(data.end(), positive.begin(), positive.end());
}  // для тестов

}  // namespace

FrolovaSRadixSortDoubleSTL::FrolovaSRadixSortDoubleSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool FrolovaSRadixSortDoubleSTL::ValidationImpl() {
  return !GetInput().empty();
}

bool FrolovaSRadixSortDoubleSTL::PreProcessingImpl() {
  return true;
}

bool FrolovaSRadixSortDoubleSTL::RunImpl() {
  const std::vector<double> &input = GetInput();
  if (input.empty()) {
    return false;
  }

  std::vector<double> working = input;
  const std::size_t n = working.size();

  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 2;
  }
  if (num_threads > n) {
    num_threads = static_cast<unsigned int>(n);
  }

  std::vector<double> temp(n);

  for (int pass = 0; pass < kNumPasses; ++pass) {
    auto histogram = ComputeHistogram(working, pass, num_threads);
    auto offsets = BuildOffsets(histogram);
    Distribute(working, temp, offsets, pass);
    working.swap(temp);
  }

  FixNegativeOrder(working);

  GetOutput() = std::move(working);
  return true;
}

bool FrolovaSRadixSortDoubleSTL::PostProcessingImpl() {
  return true;
}

}  // namespace frolova_s_radix_sort_double
