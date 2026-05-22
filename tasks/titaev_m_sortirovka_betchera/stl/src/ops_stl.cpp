#include "titaev_m_sortirovka_betchera/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <future>
#include <thread>
#include <vector>

#include "titaev_m_sortirovka_betchera/common/include/common.hpp"

namespace titaev_m_sortirovka_betchera {

namespace {

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

void RadixPass(int pass_num, size_t count, const std::vector<uint64_t> &source, std::vector<uint64_t> &dest) {
  constexpr size_t kBuckets = 256;

  std::vector<size_t> histogram(kBuckets, 0);

  for (size_t i = 0; i < count; ++i) {
    size_t bucket = (source[i] >> (static_cast<size_t>(pass_num) * 8)) & 255U;

    histogram[bucket]++;
  }

  std::vector<size_t> offsets(kBuckets, 0);

  for (size_t i = 1; i < kBuckets; ++i) {
    offsets[i] = offsets[i - 1] + histogram[i - 1];
  }

  for (size_t bucket = 0; bucket < kBuckets; ++bucket) {
    size_t pos = offsets[bucket];

    for (size_t i = 0; i < count; ++i) {
      if (((source[i] >> (static_cast<size_t>(pass_num) * 8)) & 255U) == bucket) {
        dest[pos++] = source[i];
      }
    }
  }
}

void BatcherCompareSwap(OutType &res, size_t i, size_t n, size_t step, size_t stage) {
  size_t j = i ^ stage;

  if (j <= i || j >= n) {
    return;
  }

  bool asc = (i & step) == 0;
  bool need_swap = asc ? (res[i] > res[j]) : (res[i] < res[j]);

  if (need_swap) {
    std::swap(res[i], res[j]);
  }
}

void BatcherStepRange(OutType &res, size_t begin, size_t end, size_t n, size_t step, size_t stage) {
  for (size_t i = begin; i < end; ++i) {
    BatcherCompareSwap(res, i, n, step, stage);
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
  const auto &input = GetInput();
  auto &output = GetOutput();

  const size_t size = input.size();

  output.resize(size);

  for (size_t i = 0; i < size; ++i) {
    output[i] = input[i];
  }

  return true;
}

void TitaevSortirovkaBetcheraSTL::ConvertToKeys(const InType &input, std::vector<uint64_t> &keys) {
  const size_t size = input.size();

  const size_t threads = std::max(1U, std::thread::hardware_concurrency());

  const size_t block = (size + threads - 1) / threads;

  std::vector<std::future<void>> futures;

  for (size_t ti = 0; ti < threads; ++ti) {
    const size_t begin = ti * block;

    if (begin >= size) {
      break;
    }

    const size_t end = std::min(begin + block, size);

    futures.emplace_back(std::async(std::launch::async, [&input, &keys, begin, end]() {
      for (size_t i = begin; i < end; ++i) {
        keys[i] = DoubleToOrderedUint(input[i]);
      }
    }));
  }

  for (auto &future : futures) {
    future.get();
  }
}

void TitaevSortirovkaBetcheraSTL::RadixSortParallel(std::vector<uint64_t> &keys) {
  const size_t size = keys.size();

  if (size <= 1) {
    return;
  }

  std::vector<uint64_t> temp(size);

  for (int pass = 0; pass < 8; ++pass) {
    if (pass % 2 == 0) {
      RadixPass(pass, size, keys, temp);
    } else {
      RadixPass(pass, size, temp, keys);
    }
  }
}

void TitaevSortirovkaBetcheraSTL::ConvertFromKeys(const std::vector<uint64_t> &keys, OutType &output) {
  const size_t size = keys.size();

  output.resize(size);

  const size_t threads = std::max(1U, std::thread::hardware_concurrency());

  const size_t block = (size + threads - 1) / threads;

  std::vector<std::future<void>> futures;

  for (size_t ti = 0; ti < threads; ++ti) {
    const size_t begin = ti * block;

    if (begin >= size) {
      break;
    }

    const size_t end = std::min(begin + block, size);

    futures.emplace_back(std::async(std::launch::async, [&keys, &output, begin, end]() {
      for (size_t i = begin; i < end; ++i) {
        output[i] = OrderedUintToDouble(keys[i]);
      }
    }));
  }

  for (auto &future : futures) {
    future.get();
  }
}

void TitaevSortirovkaBetcheraSTL::BatcherStepParallel(OutType &res, size_t n, size_t step, size_t stage) {
  const size_t threads = std::max(1U, std::thread::hardware_concurrency());

  const size_t block = (n + threads - 1) / threads;

  std::vector<std::future<void>> futures;

  for (size_t ti = 0; ti < threads; ++ti) {
    const size_t begin = ti * block;

    if (begin >= n) {
      break;
    }

    const size_t end = std::min(begin + block, n);

    futures.emplace_back(std::async(std::launch::async, [&res, begin, end, n, step, stage]() {
      BatcherStepRange(res, begin, end, n, step, stage);
    }));
  }

  for (auto &future : futures) {
    future.get();
  }
}

void TitaevSortirovkaBetcheraSTL::BatcherSortParallel() {
  auto &res = GetOutput();

  const size_t n = res.size();

  for (size_t step = 1; step < n; step <<= 1) {
    for (size_t stage = step; stage > 0; stage >>= 1) {
      BatcherStepParallel(res, n, step, stage);
    }
  }
}

bool TitaevSortirovkaBetcheraSTL::RunImpl() {
  const auto &input = GetInput();

  const size_t size = input.size();

  auto &output = GetOutput();

  output.resize(size);

  if (size <= 1) {
    for (size_t i = 0; i < size; ++i) {
      output[i] = input[i];
    }
    return true;
  }

  std::vector<uint64_t> keys(size);

  ConvertToKeys(input, keys);

  RadixSortParallel(keys);

  ConvertFromKeys(keys, output);

  if ((size & (size - 1)) == 0) {
    BatcherSortParallel();
  }

  return true;
}

bool TitaevSortirovkaBetcheraSTL::PostProcessingImpl() {
  return true;
}

}  // namespace titaev_m_sortirovka_betchera
