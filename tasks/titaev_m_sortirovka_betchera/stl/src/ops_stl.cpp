#include "titaev_m_sortirovka_betchera/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <thread>
#include <vector>

#include "titaev_m_sortirovka_betchera/common/include/common.hpp"

namespace titaev_m_sortirovka_betchera {

namespace {

uint64_t DoubleToOrderedUint(double value) {
  uint64_t x = 0;
  std::memcpy(&x, &value, sizeof(double));
  constexpr uint64_t kSignMask = (1ULL << 63);
  if ((x & kSignMask) != 0ULL) {
    x = ~x;
  } else {
    x ^= kSignMask;
  }
  return x;
}

double OrderedUintToDouble(uint64_t x) {
  constexpr uint64_t kSignMask = (1ULL << 63);
  if ((x & kSignMask) != 0ULL) {
    x ^= kSignMask;
  } else {
    x = ~x;
  }
  double result = 0.0;
  std::memcpy(&result, &x, sizeof(double));
  return result;
}

unsigned int GetThreadCount() {
  unsigned int hw = std::thread::hardware_concurrency();
  return hw == 0 ? 1U : hw;
}
}  // namespace

TitaevSortirovkaBetcheraSTL::TitaevSortirovkaBetcheraSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool TitaevSortirovkaBetcheraSTL::ValidationImpl() {
  return !GetInput().empty();
}

bool TitaevSortirovkaBetcheraSTL::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

void TitaevSortirovkaBetcheraSTL::ConvertToKeys(const InType &input, std::vector<uint64_t> &keys) {
  const size_t n = input.size();
  const unsigned int num_threads = GetThreadCount();
  std::vector<std::thread> threads;
  threads.reserve(num_threads);
  const size_t chunk = (n + num_threads - 1) / num_threads;

  for (unsigned int t = 0; t < num_threads; t++) {
    const size_t begin = t * chunk;
    const size_t end = std::min(begin + chunk, n);
    if (begin >= end) {
      break;
    }
    threads.emplace_back([&input, &keys, begin, end]() {
      for (size_t i = begin; i < end; i++) {
        keys[i] = DoubleToOrderedUint(input[i]);
      }
    });
  }
  for (auto &th : threads) {
    th.join();
  }
}

void TitaevSortirovkaBetcheraSTL::RadixSort(std::vector<uint64_t> &keys) {
  const size_t n = keys.size();
  if (n <= 1) {
    return;
  }

  constexpr int kBits = 8;
  constexpr int kBuckets = 1 << kBits;
  constexpr int kPasses = 64 / kBits;

  const unsigned int num_threads = GetThreadCount();
  std::vector<uint64_t> tmp(n);

  for (int pass = 0; pass < kPasses; pass++) {
    std::vector<std::vector<size_t>> local_count(num_threads, std::vector<size_t>(kBuckets, 0));
    std::vector<std::thread> threads;
    threads.reserve(num_threads);
    const size_t chunk = (n + num_threads - 1) / num_threads;

    for (unsigned int t = 0; t < num_threads; t++) {
      const size_t begin = t * chunk;
      const size_t end = std::min(begin + chunk, n);
      if (begin >= end) {
        break;
      }
      threads.emplace_back([&keys, &local_count, t, begin, end, pass]() {
        auto &lc = local_count[t];
        for (size_t i = begin; i < end; i++) {
          size_t bucket = (keys[i] >> (pass * kBits)) & (kBuckets - 1);
          lc[bucket]++;
        }
      });
    }
    for (auto &th : threads) {
      th.join();
    }

    std::vector<size_t> count(kBuckets, 0);
    for (int b = 0; b < kBuckets; b++) {
      for (unsigned int t = 0; t < num_threads; t++) {
        count[b] += local_count[t][b];
      }
    }
    for (int i = 1; i < kBuckets; i++) {
      count[i] += count[i - 1];
    }

    for (size_t i = n; i-- > 0;) {
      size_t bucket = (keys[i] >> (pass * kBits)) & (kBuckets - 1);
      tmp[--count[bucket]] = keys[i];
    }

    keys.swap(tmp);
  }
}

void TitaevSortirovkaBetcheraSTL::ConvertFromKeys(const std::vector<uint64_t> &keys, OutType &output) {
  const size_t n = keys.size();
  output.resize(n);
  const unsigned int num_threads = GetThreadCount();
  std::vector<std::thread> threads;
  threads.reserve(num_threads);
  const size_t chunk = (n + num_threads - 1) / num_threads;

  for (unsigned int t = 0; t < num_threads; t++) {
    const size_t begin = t * chunk;
    const size_t end = std::min(begin + chunk, n);
    if (begin >= end) {
      break;
    }
    threads.emplace_back([&keys, &output, begin, end]() {
      for (size_t i = begin; i < end; i++) {
        output[i] = OrderedUintToDouble(keys[i]);
      }
    });
  }
  for (auto &th : threads) {
    th.join();
  }
}

void TitaevSortirovkaBetcheraSTL::BatcherSort() {
  auto &result = GetOutput();
  const size_t n = result.size();
  if (n < 2) {
    return;
  }

  const unsigned int num_threads = GetThreadCount();

  for (size_t k = 2; k <= n; k <<= 1) {
    for (size_t j = k >> 1; j > 0; j >>= 1) {
      std::vector<std::thread> threads;
      threads.reserve(num_threads);
      const size_t chunk = (n + num_threads - 1) / num_threads;

      for (unsigned int t = 0; t < num_threads; t++) {
        const size_t begin = t * chunk;
        const size_t end = std::min(begin + chunk, n);
        if (begin >= end) {
          break;
        }
        threads.emplace_back([&result, n, k, j, begin, end]() {
          for (size_t i = begin; i < end; i++) {
            const size_t l = i ^ j;
            if (l > i && l < n) {
              const bool ascending = ((i & k) == 0);
              const bool need_swap = ascending ? (result[i] > result[l]) : (result[i] < result[l]);
              if (need_swap) {
                std::swap(result[i], result[l]);
              }
            }
          }
        });
      }
      for (auto &th : threads) {
        th.join();
      }
    }
  }
}

bool TitaevSortirovkaBetcheraSTL::RunImpl() {
  auto &input = GetInput();
  const size_t n = input.size();
  if (n <= 1) {
    return true;
  }
  std::vector<uint64_t> keys(n);
  ConvertToKeys(input, keys);
  RadixSort(keys);
  ConvertFromKeys(keys, GetOutput());
  if ((n & (n - 1)) == 0) {
    BatcherSort();
  }
  return true;
}

bool TitaevSortirovkaBetcheraSTL::PostProcessingImpl() {
  return true;
}

}  // namespace titaev_m_sortirovka_betchera
