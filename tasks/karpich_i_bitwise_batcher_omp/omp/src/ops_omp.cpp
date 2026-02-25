#include "karpich_i_bitwise_batcher_omp/omp/include/ops_omp.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <utility>
#include <vector>

#include "karpich_i_bitwise_batcher_omp/common/include/common.hpp"
#include "util/include/util.hpp"

namespace karpich_i_bitwise_batcher_omp {

namespace {

constexpr int kBytesPerInt = 4;
constexpr int kBitsPerByte = 8;
constexpr int kRadixSize = 256;
constexpr uint32_t kSignBitMask = 0x80000000U;
constexpr uint32_t kByteMask = 0xFFU;

void RadixSortRange(std::vector<int> &arr, int lo, int hi) {
  int n = hi - lo;
  if (n <= 1) {
    return;
  }

  std::vector<int> buf(n);

  for (int byte_idx = 0; byte_idx < kBytesPerInt; byte_idx++) {
    int shift = byte_idx * kBitsPerByte;
    std::array<int, kRadixSize> count{};

    for (int i = lo; i < hi; i++) {
      auto val = static_cast<uint32_t>(arr[i]);
      if (byte_idx == kBytesPerInt - 1) {
        val ^= kSignBitMask;
      }
      count.at((val >> shift) & kByteMask)++;
    }

    int prefix = 0;
    for (int &c : count) {
      int tmp = c;
      c = prefix;
      prefix += tmp;
    }

    for (int i = lo; i < hi; i++) {
      auto val = static_cast<uint32_t>(arr[i]);
      if (byte_idx == kBytesPerInt - 1) {
        val ^= kSignBitMask;
      }
      buf.at(count.at((val >> shift) & kByteMask)++) = arr[i];
    }

    std::copy(buf.begin(), buf.begin() + n, arr.begin() + lo);
  }
}

void CompareExchange(std::vector<int> &data, int pos_a, int pos_b, int len) {
  std::vector<int> merged(static_cast<size_t>(2) * len);
  std::merge(data.begin() + pos_a, data.begin() + pos_a + len, data.begin() + pos_b, data.begin() + pos_b + len,
             merged.begin());
  std::copy(merged.begin(), merged.begin() + len, data.begin() + pos_a);
  std::copy(merged.begin() + len, merged.end(), data.begin() + pos_b);
}

std::vector<std::pair<int, int>> BuildBatcherNetwork(int n) {
  std::vector<std::pair<int, int>> net;
  for (int pw = 1; pw < n; pw *= 2) {
    for (int k = pw; k >= 1; k /= 2) {
      for (int j = k % pw; j <= n - 1 - k; j += 2 * k) {
        for (int i = 0; i < std::min(k, n - j - k); i++) {
          if ((j + i) / (2 * pw) == (j + i + k) / (2 * pw)) {
            net.emplace_back(j + i, j + i + k);
          }
        }
      }
    }
  }
  return net;
}

}  // namespace

KarpichIBitwiseBatcherOMP::KarpichIBitwiseBatcherOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

bool KarpichIBitwiseBatcherOMP::ValidationImpl() {
  return GetInput() > 0 && GetOutput() == 0;
}

bool KarpichIBitwiseBatcherOMP::PreProcessingImpl() {
  int n = GetInput();
  data_.resize(n);
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int> dist(0, n);
  for (int i = 0; i < n; i++) {
    data_[i] = dist(gen);
  }
  return true;
}

bool KarpichIBitwiseBatcherOMP::RunImpl() {
  int n = static_cast<int>(data_.size());
  if (n <= 1) {
    return true;
  }

  int num_thr = ppc::util::GetNumThreads();
  if (num_thr > n) {
    num_thr = 1;
  }

  int chunk = (n + num_thr - 1) / num_thr;
  int padded = chunk * num_thr;
  data_.resize(padded, std::numeric_limits<int>::max());

  auto &data_ref = data_;
#pragma omp parallel for default(none) shared(data_ref, chunk, num_thr) num_threads(num_thr)
  for (int tid = 0; tid < num_thr; tid++) {
    RadixSortRange(data_ref, tid * chunk, (tid + 1) * chunk);
  }

  auto net = BuildBatcherNetwork(num_thr);
  for (auto &[a, b] : net) {
    CompareExchange(data_, a * chunk, b * chunk, chunk);
  }

  data_.resize(n);
  return true;
}

bool KarpichIBitwiseBatcherOMP::PostProcessingImpl() {
  GetOutput() = GetInput();
  return std::ranges::is_sorted(data_);
}

}  // namespace karpich_i_bitwise_batcher_omp
