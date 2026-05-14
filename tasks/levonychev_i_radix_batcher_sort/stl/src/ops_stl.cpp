#include "levonychev_i_radix_batcher_sort/stl/include/ops_stl.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <future>
#include <iterator>
#include <thread>
#include <vector>

#include "levonychev_i_radix_batcher_sort/common/include/common.hpp"
namespace levonychev_i_radix_batcher_sort {

LevonychevIRadixBatcherSortSTL::LevonychevIRadixBatcherSortSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

void LevonychevIRadixBatcherSortSTL::RadixSortSequential(std::vector<int> &arr) {
  if (arr.empty()) {
    return;
  }
  int n = static_cast<int>(arr.size());
  std::vector<int> buffer(n);

  for (int byte_idx = 0; byte_idx < 4; ++byte_idx) {
    std::array<int, 256> count{};
    bool is_last_byte = (byte_idx == 3);

    for (int x : arr) {
      unsigned char b = (static_cast<unsigned int>(x) >> (byte_idx * 8)) & 0xFF;
      if (is_last_byte) {
        b ^= 0x80;
      }
      count.at(b)++;
    }

    std::array<size_t, 256> offsets{};
    offsets.at(0) = 0;
    for (int i = 1; i < 256; ++i) {
      offsets.at(i) = offsets.at(i - 1) + count.at(i - 1);
    }

    for (int x : arr) {
      unsigned char b = (static_cast<unsigned int>(x) >> (byte_idx * 8)) & 0xFF;
      if (is_last_byte) {
        b ^= 0x80;
      }
      buffer.at(offsets.at(b)++) = x;
    }
    arr = buffer;
  }
}

void LevonychevIRadixBatcherSortSTL::MergeAndSplit(std::vector<int> &left_block, std::vector<int> &right_block) {
  std::vector<int> merged;
  merged.reserve(left_block.size() + right_block.size());

  std::ranges::merge(left_block, right_block, std::back_inserter(merged));

  auto mid = static_cast<std::ptrdiff_t>(left_block.size());
  std::copy(merged.begin(), merged.begin() + mid, left_block.begin());
  std::copy(merged.begin() + mid, merged.end(), right_block.begin());
}

bool LevonychevIRadixBatcherSortSTL::RunImpl() {
  const std::vector<int> data = GetInput();
  if (data.size() <= 1) {
    GetOutput() = data;
    return true;
  }

  const int num_blocks = GetNumBlocks(static_cast<int>(data.size()));
  auto blocks = DistributeData(data, num_blocks);

  ParallelRadixPhase(blocks);
  BatcherMergePhase(blocks);

  GetOutput().clear();
  GetOutput().reserve(data.size());
  for (const auto &b : blocks) {
    GetOutput().insert(GetOutput().end(), b.begin(), b.end());
  }

  return true;
}

int LevonychevIRadixBatcherSortSTL::GetNumBlocks(int n) {
  unsigned int threads_supported = std::thread::hardware_concurrency();
  int max_threads = static_cast<int>(threads_supported == 0 ? 2 : threads_supported);
  int num_blocks = 1;
  while (num_blocks * 2 <= max_threads && num_blocks * 2 <= n) {
    num_blocks *= 2;
  }
  return num_blocks;
}

std::vector<std::vector<int>> LevonychevIRadixBatcherSortSTL::DistributeData(const std::vector<int> &data,
                                                                             int num_blocks) {
  std::vector<std::vector<int>> blocks(num_blocks);
  const int n = static_cast<int>(data.size());
  const int base_size = n / num_blocks;
  const int extra = n % num_blocks;

  int current_pos = 0;
  for (int i = 0; i < num_blocks; ++i) {
    int size = base_size + (i < extra ? 1 : 0);
    blocks.at(i).assign(data.begin() + current_pos, data.begin() + current_pos + size);
    current_pos += size;
  }
  return blocks;
}

void LevonychevIRadixBatcherSortSTL::ParallelRadixPhase(std::vector<std::vector<int>> &blocks) {
  std::vector<std::future<void>> futures;
  futures.reserve(blocks.size());

  for (auto &block : blocks) {
    futures.push_back(std::async(std::launch::async, [&block]() { RadixSortSequential(block); }));
  }
  for (auto &f : futures) {
    f.wait();
  }
}

void LevonychevIRadixBatcherSortSTL::BatcherMergePhase(std::vector<std::vector<int>> &blocks) {
  const int n_blocks = static_cast<int>(blocks.size());
  for (int pk = 1; pk < n_blocks; pk <<= 1) {
    for (int k = pk; k > 0; k >>= 1) {
      BatcherMergeStep(blocks, pk, k);
    }
  }
}

void LevonychevIRadixBatcherSortSTL::BatcherMergeStep(std::vector<std::vector<int>> &blocks, int p, int k) {
  const int n_blocks = static_cast<int>(blocks.size());
  std::vector<std::future<void>> futures;
  futures.reserve(static_cast<size_t>(n_blocks));

  for (int j = k % p; j <= n_blocks - 1 - k; j += 2 * k) {
    for (int i = 0; i < std::min(k, n_blocks - j - k); ++i) {
      int idx1 = j + i;
      int idx2 = j + i + k;

      if ((idx1 / (p * 2)) == (idx2 / (p * 2))) {
        futures.push_back(std::async(std::launch::async, [&blocks, idx1, idx2]() {
          MergeAndSplit(blocks.at(static_cast<size_t>(idx1)), blocks.at(static_cast<size_t>(idx2)));
        }));
      }
    }
  }

  for (auto &f : futures) {
    f.wait();
  }
}

bool LevonychevIRadixBatcherSortSTL::ValidationImpl() {
  return !GetInput().empty();
}
bool LevonychevIRadixBatcherSortSTL::PreProcessingImpl() {
  return true;
}
bool LevonychevIRadixBatcherSortSTL::PostProcessingImpl() {
  return true;
}

}  // namespace levonychev_i_radix_batcher_sort
