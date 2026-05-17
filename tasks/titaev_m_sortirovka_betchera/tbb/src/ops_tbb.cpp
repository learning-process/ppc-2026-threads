#include "titaev_m_sortirovka_betchera/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <vector>

#include "titaev_m_sortirovka_betchera/common/include/common.hpp"

namespace titaev_m_sortirovka_betchera {

namespace {
uint64_t DoubleToOrderedUint(double value) {
  uint64_t x_val = 0;
  std::memcpy(&x_val, &value, sizeof(double));
  constexpr uint64_t kSignMask = (1ULL << 63);
  if ((x_val & kSignMask) != 0ULL) {
    x_val = ~x_val;
  } else {
    x_val ^= kSignMask;
  }
  return x_val;
}

double OrderedUintToDouble(uint64_t x_val) {
  constexpr uint64_t kSignMask = (1ULL << 63);
  if ((x_val & kSignMask) != 0ULL) {
    x_val ^= kSignMask;
  } else {
    x_val = ~x_val;
  }
  double res_val = 0.0;
  std::memcpy(&res_val, &x_val, sizeof(double));
  return res_val;
}

void RadixPass(int pass_num, size_t count_n, const std::vector<uint64_t> &source, std::vector<uint64_t> &dest) {
  constexpr size_t kNumBuckets = 256;
  using LocalHistogram = std::vector<size_t>;
  tbb::enumerable_thread_specific<LocalHistogram> histograms(LocalHistogram(kNumBuckets, 0));

  tbb::parallel_for(tbb::blocked_range<size_t>(0, count_n), [&](const tbb::blocked_range<size_t> &range) {
    auto &my_hist = histograms.local();
    for (size_t i = range.begin(); i < range.end(); ++i) {
      size_t bucket_idx = (source[i] >> (static_cast<size_t>(pass_num) * 8)) & 255;
      my_hist[bucket_idx]++;
    }
  });

  std::vector<size_t> global_hist(kNumBuckets, 0);
  for (const auto &hist : histograms) {
    for (size_t b_idx = 0; b_idx < kNumBuckets; ++b_idx) {
      global_hist[b_idx] += hist[b_idx];
    }
  }

  std::vector<size_t> start_offsets(kNumBuckets, 0);
  for (size_t b_idx = 1; b_idx < kNumBuckets; ++b_idx) {
    start_offsets[b_idx] = start_offsets[b_idx - 1] + global_hist[b_idx - 1];
  }

  for (size_t b_idx = 0; b_idx < kNumBuckets; ++b_idx) {
    size_t write_pos = start_offsets[b_idx];
    for (size_t i = 0; i < count_n; ++i) {
      if (((source[i] >> (static_cast<size_t>(pass_num) * 8)) & 255) == b_idx) {
        dest[write_pos++] = source[i];
      }
    }
  }
}
}  // namespace

TitaevSortirovkaBetcheraTBB::TitaevSortirovkaBetcheraTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool TitaevSortirovkaBetcheraTBB::ValidationImpl() {
  return !GetInput().empty();
}
bool TitaevSortirovkaBetcheraTBB::PreProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

void TitaevSortirovkaBetcheraTBB::ConvertToKeys(const InType &input, std::vector<uint64_t> &keys) {
  tbb::parallel_for(tbb::blocked_range<size_t>(0, input.size()), [&](const tbb::blocked_range<size_t> &range) {
    for (size_t i = range.begin(); i < range.end(); ++i) {
      keys[i] = DoubleToOrderedUint(input[i]);
    }
  });
}

void TitaevSortirovkaBetcheraTBB::RadixSortParallel(std::vector<uint64_t> &keys) {
  const size_t count_n = keys.size();
  if (count_n <= 1) {
    return;
  }
  std::vector<uint64_t> temp_keys(count_n);
  for (int pass_idx = 0; pass_idx < 8; ++pass_idx) {
    if (pass_idx % 2 == 0) {
      RadixPass(pass_idx, count_n, keys, temp_keys);
    } else {
      RadixPass(pass_idx, count_n, temp_keys, keys);
    }
  }
}

void TitaevSortirovkaBetcheraTBB::ConvertFromKeys(const std::vector<uint64_t> &keys, OutType &output) {
  output.resize(keys.size());
  tbb::parallel_for(tbb::blocked_range<size_t>(0, keys.size()), [&](const tbb::blocked_range<size_t> &range) {
    for (size_t i = range.begin(); i < range.end(); ++i) {
      output[i] = OrderedUintToDouble(keys[i]);
    }
  });
}

void TitaevSortirovkaBetcheraTBB::BatcherStepParallel(OutType &res_vec, size_t total_n, size_t step_size,
                                                      size_t stage_dist) {
  tbb::parallel_for(tbb::blocked_range<size_t>(0, total_n), [&](const tbb::blocked_range<size_t> &range) {
    for (size_t i = range.begin(); i < range.end(); ++i) {
      size_t j_idx = i ^ stage_dist;
      if (j_idx > i && j_idx < total_n) {
        bool is_asc = (i & step_size) == 0;
        if (is_asc ? (res_vec[i] > res_vec[j_idx]) : (res_vec[i] < res_vec[j_idx])) {
          std::swap(res_vec[i], res_vec[j_idx]);
        }
      }
    }
  });
}

void TitaevSortirovkaBetcheraTBB::BatcherSortParallel() {
  auto &res_vec = GetOutput();
  const size_t total_n = res_vec.size();
  for (size_t step_size = 1; step_size < total_n; step_size <<= 1) {
    for (size_t stage_dist = step_size; stage_dist > 0; stage_dist >>= 1) {
      BatcherStepParallel(res_vec, total_n, step_size, stage_dist);
    }
  }
}

bool TitaevSortirovkaBetcheraTBB::RunImpl() {
  auto &input_data = GetInput();
  const size_t data_size = input_data.size();
  if (data_size <= 1) {
    return true;
  }

  std::vector<uint64_t> keys_vec(data_size);
  ConvertToKeys(input_data, keys_vec);
  RadixSortParallel(keys_vec);
  ConvertFromKeys(keys_vec, GetOutput());

  if ((data_size & (data_size - 1)) == 0) {
    BatcherSortParallel();
  }
  return true;
}

bool TitaevSortirovkaBetcheraTBB::PostProcessingImpl() {
  return true;
}

}  // namespace titaev_m_sortirovka_betchera
