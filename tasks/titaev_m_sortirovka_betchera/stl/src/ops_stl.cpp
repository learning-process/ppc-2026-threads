#include "titaev_m_sortirovka_betchera/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <functional>
#include <thread>
#include <vector>

#include "titaev_m_sortirovka_betchera/common/include/common.hpp"

namespace titaev_m_sortirovka_betchera {

namespace {

constexpr int kBits = 8;
constexpr int kBuckets = 1 << kBits;
constexpr int kPasses = 64 / kBits;
constexpr std::size_t kSequentialThreshold = 1 << 16;

uint64_t DoubleToOrderedUint(double value) {
  uint64_t bits = 0;
  std::memcpy(&bits, &value, sizeof(double));
  constexpr uint64_t kSignMask = (1ULL << 63);
  if ((bits & kSignMask) != 0ULL) {
    bits = ~bits;
  } else {
    bits ^= kSignMask;
  }
  return bits;
}

double OrderedUintToDouble(uint64_t bits) {
  constexpr uint64_t kSignMask = (1ULL << 63);
  if ((bits & kSignMask) != 0ULL) {
    bits ^= kSignMask;
  } else {
    bits = ~bits;
  }
  double result = 0.0;
  std::memcpy(&result, &bits, sizeof(double));
  return result;
}

unsigned int GetThreadCount() {
  const unsigned int hw = std::thread::hardware_concurrency();
  return hw == 0 ? 1U : hw;
}

void RunInParallel(std::size_t total, const std::function<void(std::size_t, std::size_t)> &body) {
  const unsigned int num_threads = (total < kSequentialThreshold) ? 1U : GetThreadCount();

  if (num_threads <= 1) {
    body(0, total);
    return;
  }

  std::vector<std::thread> threads;
  threads.reserve(num_threads);
  const std::size_t chunk = (total + num_threads - 1) / num_threads;

  for (unsigned int thr = 0; thr < num_threads; thr++) {
    const std::size_t begin = thr * chunk;
    const std::size_t end = std::min(begin + chunk, total);
    if (begin >= end) {
      break;
    }
    threads.emplace_back([&body, begin, end]() { body(begin, end); });
  }
  for (auto &th : threads) {
    th.join();
  }
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
  RunInParallel(input.size(), [&input, &keys](std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      keys[i] = DoubleToOrderedUint(input[i]);
    }
  });
}

void TitaevSortirovkaBetcheraSTL::CountSequential(const std::vector<uint64_t> &keys, std::vector<std::size_t> &count,
                                                  int pass) {
  const std::size_t n = keys.size();
  for (std::size_t i = 0; i < n; i++) {
    const std::size_t bucket = (keys[i] >> (pass * kBits)) & (kBuckets - 1);
    count[bucket]++;
  }
}

void TitaevSortirovkaBetcheraSTL::CountParallel(const std::vector<uint64_t> &keys, std::vector<std::size_t> &count,
                                                int pass, unsigned int num_threads) {
  const std::size_t n = keys.size();
  std::vector<std::vector<std::size_t>> local_count(num_threads, std::vector<std::size_t>(kBuckets, 0));
  const std::size_t chunk = (n + num_threads - 1) / num_threads;

  std::vector<std::thread> threads;
  threads.reserve(num_threads);
  for (unsigned int thr = 0; thr < num_threads; thr++) {
    const std::size_t begin = thr * chunk;
    const std::size_t end = std::min(begin + chunk, n);
    if (begin >= end) {
      break;
    }
    threads.emplace_back([&keys, &local_count, thr, begin, end, pass]() {
      auto &lc = local_count[thr];
      for (std::size_t i = begin; i < end; i++) {
        const std::size_t bucket = (keys[i] >> (pass * kBits)) & (kBuckets - 1);
        lc[bucket]++;
      }
    });
  }
  for (auto &th : threads) {
    th.join();
  }

  for (int bucket_idx = 0; bucket_idx < kBuckets; bucket_idx++) {
    for (unsigned int thr = 0; thr < num_threads; thr++) {
      count[bucket_idx] += local_count[thr][bucket_idx];
    }
  }
}

void TitaevSortirovkaBetcheraSTL::RadixCountPass(std::vector<uint64_t> &keys, std::vector<uint64_t> &tmp, int pass) {
  const std::size_t n = keys.size();
  const unsigned int num_threads = (n < kSequentialThreshold) ? 1U : GetThreadCount();

  std::vector<std::size_t> count(kBuckets, 0);
  if (num_threads <= 1) {
    CountSequential(keys, count, pass);
  } else {
    CountParallel(keys, count, pass, num_threads);
  }

  for (int i = 1; i < kBuckets; i++) {
    count[i] += count[i - 1];
  }
  for (std::size_t i = n; i-- > 0;) {
    const std::size_t bucket = (keys[i] >> (pass * kBits)) & (kBuckets - 1);
    tmp[--count[bucket]] = keys[i];
  }
  keys.swap(tmp);
}
void TitaevSortirovkaBetcheraSTL::RadixSort(std::vector<uint64_t> &keys) {
  const std::size_t n = keys.size();
  if (n <= 1) {
    return;
  }
  std::vector<uint64_t> tmp(n);
  for (int pass = 0; pass < kPasses; pass++) {
    RadixCountPass(keys, tmp, pass);
  }
}

void TitaevSortirovkaBetcheraSTL::ConvertFromKeys(const std::vector<uint64_t> &keys, OutType &output) {
  output.resize(keys.size());
  RunInParallel(keys.size(), [&keys, &output](std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      output[i] = OrderedUintToDouble(keys[i]);
    }
  });
}

void TitaevSortirovkaBetcheraSTL::BatcherStage(OutType &result, std::size_t array_size, std::size_t block,
                                               std::size_t step) {
  RunInParallel(array_size, [&result, block, step](std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      const std::size_t partner = i ^ step;
      if (partner <= i) {
        continue;
      }
      const bool ascending = ((i & block) == 0);
      const bool need_swap = ascending ? (result[i] > result[partner]) : (result[i] < result[partner]);
      if (need_swap) {
        std::swap(result[i], result[partner]);
      }
    }
  });
}

void TitaevSortirovkaBetcheraSTL::BatcherSort() {
  auto &result = GetOutput();
  const std::size_t n = result.size();
  if (n < 2) {
    return;
  }
  for (std::size_t block = 2; block <= n; block <<= 1) {
    for (std::size_t step = block >> 1; step > 0; step >>= 1) {
      BatcherStage(result, n, block, step);
    }
  }
}

bool TitaevSortirovkaBetcheraSTL::RunImpl() {
  auto &input = GetInput();
  const std::size_t n = input.size();
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
