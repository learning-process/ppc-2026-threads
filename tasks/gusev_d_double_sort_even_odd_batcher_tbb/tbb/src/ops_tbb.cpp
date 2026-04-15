#include "gusev_d_double_sort_even_odd_batcher_tbb/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range.h>
#include <tbb/global_control.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_invoke.h>
#include <tbb/task_arena.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <utility>
#include <vector>

#include "gusev_d_double_sort_even_odd_batcher_tbb/common/include/common.hpp"
#include "util/include/util.hpp"

namespace gusev_d_double_sort_even_odd_batcher_tbb_task_threads {
namespace {

constexpr int kRadixPasses = 8;
constexpr int kBitsPerByte = 8;
constexpr size_t kRadixBuckets = 256;
constexpr uint64_t kBucketMask = 0xFFULL;
constexpr size_t kMinParallelElements = 128;
constexpr size_t kMinParallelPairs = 32;

static_assert((kRadixPasses % 2) == 0, "Radix sort expects the final data to remain in the input buffer");

using Block = std::vector<ValueType>;
using BlockList = std::vector<Block>;

struct BlockRange {
  size_t begin = 0;
  size_t end = 0;
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
  auto *src = &data;
  auto *dst = &buffer;

  for (int byte = 0; byte < kRadixPasses; ++byte) {
    std::array<size_t, kRadixBuckets> count{};
    const auto shift = byte * kBitsPerByte;

    for (ValueType value : *src) {
      count.at(GetBucketIndex(value, shift))++;
    }
    BuildPrefixSums(count);

    for (ValueType value : *src) {
      const auto bucket = GetBucketIndex(value, shift);
      (*dst)[count.at(bucket)++] = value;
    }

    std::swap(src, dst);
  }
}

void SplitByGlobalParity(const OutType &source, size_t global_offset, OutType &even_values, OutType &odd_values) {
  even_values.clear();
  odd_values.clear();
  even_values.reserve((source.size() + 1) / 2);
  odd_values.reserve(source.size() / 2);

  for (size_t i = 0; i < source.size(); ++i) {
    if (((global_offset + i) & 1U) == 0U) {
      even_values.push_back(source[i]);
    } else {
      odd_values.push_back(source[i]);
    }
  }
}

OutType InterleaveParityGroups(size_t total_size, const OutType &even_values, const OutType &odd_values) {
  OutType result(total_size);
  if (total_size < kMinParallelElements || even_values.empty() || odd_values.empty()) {
    for (size_t even_index = 0; even_index < even_values.size(); ++even_index) {
      result[even_index * 2] = even_values[even_index];
    }
    for (size_t odd_index = 0; odd_index < odd_values.size(); ++odd_index) {
      result[(odd_index * 2) + 1] = odd_values[odd_index];
    }
    return result;
  }

  tbb::parallel_invoke([&]() {
    for (size_t even_index = 0; even_index < even_values.size(); ++even_index) {
      result[even_index * 2] = even_values[even_index];
    }
  }, [&]() {
    for (size_t odd_index = 0; odd_index < odd_values.size(); ++odd_index) {
      result[(odd_index * 2) + 1] = odd_values[odd_index];
    }
  });

  return result;
}

void RunOddEvenPhase(OutType &data, size_t start) {
  const auto pair_count = (data.size() - start) / 2;
  if (pair_count < kMinParallelPairs) {
    for (size_t pair_index = 0; pair_index < pair_count; ++pair_index) {
      const auto left_index = start + (pair_index * 2);
      if (data[left_index] > data[left_index + 1]) {
        std::swap(data[left_index], data[left_index + 1]);
      }
    }
    return;
  }

  tbb::parallel_for(tbb::blocked_range<size_t>(0, pair_count), [&](const tbb::blocked_range<size_t> &range) {
    for (size_t pair_index = range.begin(); pair_index != range.end(); ++pair_index) {
      const auto left_index = start + (pair_index * 2);
      if (data[left_index] > data[left_index + 1]) {
        std::swap(data[left_index], data[left_index + 1]);
      }
    }
  });
}

void OddEvenFinalize(OutType &data) {
  for (size_t phase = 0; phase < data.size(); ++phase) {
    const auto start = phase & 1U;
    if (data.size() > start + 1) {
      RunOddEvenPhase(data, start);
    }
  }
}

void MergeParityGroups(const OutType &left_even, const OutType &right_even, const OutType &left_odd,
                       const OutType &right_odd, OutType &merged_even, OutType &merged_odd) {
  merged_even.clear();
  merged_odd.clear();
  merged_even.reserve(left_even.size() + right_even.size());
  merged_odd.reserve(left_odd.size() + right_odd.size());

  if ((merged_even.capacity() + merged_odd.capacity()) < kMinParallelElements) {
    std::ranges::merge(left_even, right_even, std::back_inserter(merged_even));
    std::ranges::merge(left_odd, right_odd, std::back_inserter(merged_odd));
    return;
  }

  tbb::parallel_invoke([&]() { std::ranges::merge(left_even, right_even, std::back_inserter(merged_even)); },
                       [&]() { std::ranges::merge(left_odd, right_odd, std::back_inserter(merged_odd)); });
}

void SplitBlocksByParity(const OutType &left, const OutType &right, OutType &left_even, OutType &left_odd,
                         OutType &right_even, OutType &right_odd) {
  if ((left.size() + right.size()) < kMinParallelElements) {
    SplitByGlobalParity(left, 0, left_even, left_odd);
    SplitByGlobalParity(right, left.size(), right_even, right_odd);
    return;
  }

  tbb::parallel_invoke([&]() { SplitByGlobalParity(left, 0, left_even, left_odd); },
                       [&]() { SplitByGlobalParity(right, left.size(), right_even, right_odd); });
}

OutType MergeBatcherEvenOdd(const OutType &left, const OutType &right) {
  OutType left_even;
  OutType left_odd;
  OutType right_even;
  OutType right_odd;

  SplitBlocksByParity(left, right, left_even, left_odd, right_even, right_odd);

  OutType merged_even;
  OutType merged_odd;
  MergeParityGroups(left_even, right_even, left_odd, right_odd, merged_even, merged_odd);

  auto result = InterleaveParityGroups(left.size() + right.size(), merged_even, merged_odd);
  OddEvenFinalize(result);
  return result;
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

BlockList MakeSortedBlocks(const InType &input, size_t parallelism) {
  const auto block_count = GetBlockCount(input.size(), parallelism);
  const auto total_size = input.size();

  BlockList blocks(block_count);
  if (block_count == 1 || input.size() < kMinParallelElements) {
    for (size_t block = 0; block < block_count; ++block) {
      FillAndSortBlock(input, blocks[block], GetBlockRange(block, block_count, total_size));
    }
    return blocks;
  }

  tbb::parallel_for(tbb::blocked_range<size_t>(0, block_count), [&](const tbb::blocked_range<size_t> &range) {
    for (size_t block = range.begin(); block != range.end(); ++block) {
      FillAndSortBlock(input, blocks[block], GetBlockRange(block, block_count, total_size));
    }
  });

  return blocks;
}

void MergeBlockPair(const BlockList &blocks, BlockList &next, size_t pair_index) {
  next[pair_index] = MergeBatcherEvenOdd(blocks[pair_index * 2], blocks[(pair_index * 2) + 1]);
}

BlockList MergeBlockPairs(const BlockList &blocks) {
  const auto pair_count = blocks.size() / 2;
  BlockList next((blocks.size() + 1) / 2);
  if (pair_count < kMinParallelPairs) {
    for (size_t pair = 0; pair < pair_count; ++pair) {
      MergeBlockPair(blocks, next, pair);
    }
    if ((blocks.size() & 1U) != 0U) {
      next.back() = blocks.back();
    }
    return next;
  }

  tbb::parallel_for(tbb::blocked_range<size_t>(0, pair_count), [&](const tbb::blocked_range<size_t> &range) {
    for (size_t pair = range.begin(); pair != range.end(); ++pair) {
      MergeBlockPair(blocks, next, pair);
    }
  });

  if ((blocks.size() & 1U) != 0U) {
    next.back() = blocks.back();
  }

  return next;
}

Block MergeBlocks(BlockList blocks) {
  while (blocks.size() > 1) {
    blocks = MergeBlockPairs(blocks);
  }

  return std::move(blocks.front());
}

}  // namespace

DoubleSortEvenOddBatcherTBB::DoubleSortEvenOddBatcherTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool DoubleSortEvenOddBatcherTBB::ValidationImpl() {
  return GetOutput().empty();
}

bool DoubleSortEvenOddBatcherTBB::PreProcessingImpl() {
  input_data_ = GetInput();
  result_data_.clear();
  return true;
}

bool DoubleSortEvenOddBatcherTBB::RunImpl() {
  if (input_data_.empty()) {
    result_data_.clear();
    return true;
  }

  const auto parallelism = static_cast<size_t>(std::max(1, ppc::util::GetNumThreads()));
  const tbb::global_control control(tbb::global_control::max_allowed_parallelism, parallelism);

  auto blocks = MakeSortedBlocks(input_data_, parallelism);
  result_data_ = MergeBlocks(std::move(blocks));
  return true;
}

bool DoubleSortEvenOddBatcherTBB::PostProcessingImpl() {
  GetOutput() = result_data_;
  return true;
}

}  // namespace gusev_d_double_sort_even_odd_batcher_tbb_task_threads
