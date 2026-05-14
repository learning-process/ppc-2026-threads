#include "gusev_d_double_sort_even_odd_batcher_stl/stl/include/ops_stl.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <thread>
#include <utility>
#include <vector>

#include "gusev_d_double_sort_even_odd_batcher_stl/common/include/common.hpp"
#include "util/include/util.hpp"

namespace gusev_d_double_sort_even_odd_batcher_stl_task_threads {
namespace {

constexpr int kRadixPasses = 8;
constexpr int kBitsPerByte = 8;
constexpr size_t kRadixBuckets = 256;
constexpr uint64_t kBucketMask = 0xFFULL;
constexpr size_t kMinParallelElements = 128;
constexpr size_t kMinThreadedTasks = 2;

static_assert((kRadixPasses % 2) == 0, "Radix sort expects the final data to remain in the input buffer");

using Block = std::vector<ValueType>;
using BlockList = std::vector<Block>;

struct BlockRange {
  size_t begin = 0;
  size_t end = 0;
};

struct MergeItem {
  ValueType value = 0.0;
  bool is_padding = true;
};

uint64_t DoubleToSortableKey(ValueType value) {
  const auto bits = std::bit_cast<uint64_t>(value);
  const auto sign_mask = uint64_t{1} << 63;
  return (bits & sign_mask) == 0 ? bits ^ sign_mask : ~bits;
}

size_t GetBucketIndex(ValueType value, int shift) {
  return static_cast<size_t>((DoubleToSortableKey(value) >> shift) & kBucketMask);
}

void BuildPrefixSums(std::array<size_t, kRadixBuckets> &count) {
  size_t prefix = 0;
  for (auto &value : count) {
    const auto current = value;
    value = prefix;
    prefix += current;
  }
}

void RadixSortDoubles(OutType &data) {
  if (data.size() < 2) {
    return;
  }

  OutType buffer(data.size());
  auto *source = &data;
  auto *destination = &buffer;

  for (int byte = 0; byte < kRadixPasses; ++byte) {
    std::array<size_t, kRadixBuckets> count{};
    const auto shift = byte * kBitsPerByte;

    for (ValueType value : *source) {
      count[GetBucketIndex(value, shift)]++;
    }
    BuildPrefixSums(count);

    for (ValueType value : *source) {
      const auto bucket = GetBucketIndex(value, shift);
      (*destination)[count[bucket]++] = value;
    }

    std::swap(source, destination);
  }
}

size_t NextPowerOfTwo(size_t value) {
  size_t result = 1;
  while (result < value) {
    result <<= 1U;
  }
  return result;
}

bool IsGreater(const MergeItem &lhs, const MergeItem &rhs) {
  if (lhs.is_padding != rhs.is_padding) {
    return lhs.is_padding;
  }
  return !lhs.is_padding && lhs.value > rhs.value;
}

void CompareExchange(std::vector<MergeItem> &data, size_t left, size_t right) {
  if (IsGreater(data[left], data[right])) {
    std::swap(data[left], data[right]);
  }
}

void OddEvenMerge(std::vector<MergeItem> &data, size_t first, size_t length, size_t distance) {
  const auto next_distance = distance * 2;
  if (next_distance < length) {
    OddEvenMerge(data, first, length, next_distance);
    OddEvenMerge(data, first + distance, length, next_distance);
    for (size_t i = first + distance; (i + distance) < (first + length); i += next_distance) {
      CompareExchange(data, i, i + distance);
    }
  } else {
    CompareExchange(data, first, first + distance);
  }
}

void CopyBlockToMergeItems(const Block &block, std::vector<MergeItem> &items, size_t offset) {
  for (size_t i = 0; i < block.size(); ++i) {
    items[offset + i] = {.value = block[i], .is_padding = false};
  }
}

Block ExtractMergedValues(const std::vector<MergeItem> &items, size_t result_size) {
  Block result;
  result.reserve(result_size);
  for (const auto &item : items) {
    if (!item.is_padding) {
      result.push_back(item.value);
    }
  }
  return result;
}

Block MergeBatcherEvenOdd(const Block &left, const Block &right) {
  if (left.empty()) {
    return right;
  }
  if (right.empty()) {
    return left;
  }

  const auto half_size = NextPowerOfTwo(std::max(left.size(), right.size()));
  std::vector<MergeItem> items(half_size * 2);

  CopyBlockToMergeItems(left, items, 0);
  CopyBlockToMergeItems(right, items, half_size);
  OddEvenMerge(items, 0, items.size(), 1);
  return ExtractMergedValues(items, left.size() + right.size());
}

BlockRange GetBlockRange(size_t block_index, size_t block_count, size_t total_size) {
  return {
      .begin = (block_index * total_size) / block_count,
      .end = ((block_index + 1) * total_size) / block_count,
  };
}

size_t GetBlockCount(size_t input_size, size_t parallelism) {
  return std::max<size_t>(1, std::min(input_size, parallelism));
}

void FillAndSortBlock(const InType &input, Block &block, BlockRange range) {
  block.assign(input.begin() + static_cast<std::ptrdiff_t>(range.begin),
               input.begin() + static_cast<std::ptrdiff_t>(range.end));
  RadixSortDoubles(block);
}

template <typename Function>
void RunThreadedByIndex(size_t work_count, size_t max_threads, Function &&function) {
  if (work_count < kMinThreadedTasks || max_threads < kMinThreadedTasks) {
    for (size_t index = 0; index < work_count; ++index) {
      function(index);
    }
    return;
  }

  const auto thread_count = std::min(work_count, max_threads);
  std::vector<std::thread> workers;
  workers.reserve(thread_count);
  for (size_t thread_index = 0; thread_index < thread_count; ++thread_index) {
    workers.emplace_back([&, thread_index] {
      for (size_t index = thread_index; index < work_count; index += thread_count) {
        function(index);
      }
    });
  }

  for (auto &worker : workers) {
    worker.join();
  }
}

BlockList MakeSortedBlocks(const InType &input, size_t parallelism) {
  const auto block_count = GetBlockCount(input.size(), parallelism);
  const auto total_size = input.size();
  BlockList blocks(block_count);

  const auto sort_block = [&](size_t block_index) {
    FillAndSortBlock(input, blocks[block_index], GetBlockRange(block_index, block_count, total_size));
  };

  if (input.size() < kMinParallelElements) {
    RunThreadedByIndex(block_count, 1, sort_block);
  } else {
    RunThreadedByIndex(block_count, parallelism, sort_block);
  }

  return blocks;
}

void MergeBlockPair(const BlockList &blocks, BlockList &next, size_t pair_index) {
  next[pair_index] = MergeBatcherEvenOdd(blocks[pair_index * 2], blocks[(pair_index * 2) + 1]);
}

BlockList MergeBlockPairs(const BlockList &blocks, size_t parallelism) {
  const auto pair_count = blocks.size() / 2;
  BlockList next((blocks.size() + 1) / 2);

  RunThreadedByIndex(pair_count, parallelism, [&](size_t pair_index) { MergeBlockPair(blocks, next, pair_index); });

  if ((blocks.size() & 1U) != 0U) {
    next.back() = blocks.back();
  }

  return next;
}

Block MergeBlocks(BlockList blocks, size_t parallelism) {
  while (blocks.size() > 1) {
    blocks = MergeBlockPairs(blocks, parallelism);
  }

  return std::move(blocks.front());
}

}  // namespace

DoubleSortEvenOddBatcherSTL::DoubleSortEvenOddBatcherSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool DoubleSortEvenOddBatcherSTL::ValidationImpl() {
  return GetOutput().empty();
}

bool DoubleSortEvenOddBatcherSTL::PreProcessingImpl() {
  input_data_ = GetInput();
  result_data_.clear();
  return true;
}

bool DoubleSortEvenOddBatcherSTL::RunImpl() {
  if (input_data_.empty()) {
    result_data_.clear();
    return true;
  }

  const auto parallelism = static_cast<size_t>(std::max(1, ppc::util::GetNumThreads()));
  auto blocks = MakeSortedBlocks(input_data_, parallelism);
  result_data_ = MergeBlocks(std::move(blocks), parallelism);
  return true;
}

bool DoubleSortEvenOddBatcherSTL::PostProcessingImpl() {
  GetOutput() = result_data_;
  return true;
}

}  // namespace gusev_d_double_sort_even_odd_batcher_stl_task_threads
