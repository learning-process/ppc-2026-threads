#include "gusev_d_double_sort_even_odd_batcher_omp/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <iterator>
#include <utility>
#include <vector>

#include "gusev_d_double_sort_even_odd_batcher_omp/common/include/common.hpp"

namespace gusev_d_double_sort_even_odd_batcher_omp_task_threads {
namespace {

constexpr int kRadixPasses = 8;
constexpr int kBitsPerByte = 8;
constexpr size_t kRadixBuckets = 256;
constexpr uint64_t kBucketMask = 0xFFULL;

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

void RadixSortDoubles(Block &data) {
  if (data.size() < 2) {
    return;
  }

  Block buffer(data.size());
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

void SplitByGlobalParity(const Block &source, size_t global_offset, Block &even, Block &odd) {
  even.clear();
  odd.clear();
  even.reserve((source.size() + 1) / 2);
  odd.reserve(source.size() / 2);

  for (size_t i = 0; i < source.size(); ++i) {
    if (((global_offset + i) & 1U) == 0U) {
      even.push_back(source[i]);
    } else {
      odd.push_back(source[i]);
    }
  }
}

Block InterleaveParityGroups(size_t total_size, const Block &even, const Block &odd) {
  Block result(total_size);
  size_t even_index = 0;
  size_t odd_index = 0;

  for (size_t i = 0; i < total_size; ++i) {
    if ((i & 1U) == 0U) {
      result[i] = even[even_index++];
    } else {
      result[i] = odd[odd_index++];
    }
  }

  return result;
}

void OddEvenFinalize(Block &result) {
  for (size_t phase = 0; phase < result.size(); ++phase) {
    const auto start = phase & 1U;
    for (size_t i = start; i + 1 < result.size(); i += 2) {
      if (result[i] > result[i + 1]) {
        std::swap(result[i], result[i + 1]);
      }
    }
  }
}

void SplitBlocksByParity(const Block &left, const Block &right, Block &left_even, Block &left_odd, Block &right_even,
                         Block &right_odd) {
  SplitByGlobalParity(left, 0, left_even, left_odd);
  SplitByGlobalParity(right, left.size(), right_even, right_odd);
}

void MergeParityGroups(const Block &left_even, const Block &right_even, const Block &left_odd, const Block &right_odd,
                       Block &merged_even, Block &merged_odd) {
  merged_even.clear();
  merged_odd.clear();
  merged_even.reserve(left_even.size() + right_even.size());
  merged_odd.reserve(left_odd.size() + right_odd.size());

  std::ranges::merge(left_even, right_even, std::back_inserter(merged_even));
  std::ranges::merge(left_odd, right_odd, std::back_inserter(merged_odd));
}

Block MergeBatcherEvenOdd(const Block &left, const Block &right) {
  Block left_even;
  Block left_odd;
  Block right_even;
  Block right_odd;

  SplitBlocksByParity(left, right, left_even, left_odd, right_even, right_odd);

  Block merged_even;
  Block merged_odd;
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

size_t GetBlockCount(size_t input_size) {
  const auto omp_threads = static_cast<size_t>(std::max(1, omp_get_max_threads()));
  return std::max<size_t>(1, std::min(input_size, omp_threads));
}

void FillAndSortBlock(const std::vector<ValueType> &input, Block &block, BlockRange range) {
  block.assign(input.begin() + static_cast<std::ptrdiff_t>(range.begin),
               input.begin() + static_cast<std::ptrdiff_t>(range.end));
  RadixSortDoubles(block);
}

void StoreFirstException(std::atomic_bool &has_exception, std::exception_ptr &exception) {
#pragma omp critical
  if (!has_exception.exchange(true)) {
    exception = std::current_exception();
  }
}

template <class Func>
void ExecuteWithExceptionCapture(std::atomic_bool &has_exception, std::exception_ptr &exception, Func &&func) {
  if (has_exception.load()) {
    return;
  }

  try {
    std::forward<Func>(func)();
  } catch (...) {
    StoreFirstException(has_exception, exception);
  }
}

void RethrowIfCaptured(const std::exception_ptr &exception) {
  if (exception != nullptr) {
    std::rethrow_exception(exception);
  }
}

BlockList MakeSortedBlocks(const std::vector<ValueType> &input) {
  const auto block_count = GetBlockCount(input.size());
  const auto total_size = input.size();
  const auto signed_block_count = static_cast<std::ptrdiff_t>(block_count);

  BlockList blocks(block_count);
  std::exception_ptr exception;
  std::atomic_bool has_exception{false};

#pragma omp parallel for default(none) shared(input, blocks, exception, has_exception) \
    firstprivate(block_count, signed_block_count, total_size) schedule(static) if (block_count > 1)
  for (std::ptrdiff_t block = 0; block < signed_block_count; ++block) {
    ExecuteWithExceptionCapture(has_exception, exception, [&]() {
      const auto index = static_cast<size_t>(block);
      FillAndSortBlock(input, blocks[index], GetBlockRange(index, block_count, total_size));
    });
  }

  RethrowIfCaptured(exception);

  return blocks;
}

void MergeBlockPair(const BlockList &blocks, BlockList &next, size_t pair_index) {
  next[pair_index] = MergeBatcherEvenOdd(blocks[pair_index * 2], blocks[(pair_index * 2) + 1]);
}

BlockList MergeBlockPairs(const BlockList &blocks) {
  const auto signed_pair_count = static_cast<std::ptrdiff_t>(blocks.size() / 2);
  BlockList next((blocks.size() + 1) / 2);
  std::exception_ptr exception;
  std::atomic_bool has_exception{false};

#pragma omp parallel for default(none) shared(blocks, next, exception, has_exception) firstprivate(signed_pair_count) \
    schedule(static) if (signed_pair_count > 1)
  for (std::ptrdiff_t pair = 0; pair < signed_pair_count; ++pair) {
    ExecuteWithExceptionCapture(has_exception, exception,
                                [&]() { MergeBlockPair(blocks, next, static_cast<size_t>(pair)); });
  }

  RethrowIfCaptured(exception);

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

DoubleSortEvenOddBatcherOMP::DoubleSortEvenOddBatcherOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool DoubleSortEvenOddBatcherOMP::ValidationImpl() {
  return GetOutput().empty();
}

bool DoubleSortEvenOddBatcherOMP::PreProcessingImpl() {
  input_data_ = GetInput();
  result_data_.clear();
  return true;
}

bool DoubleSortEvenOddBatcherOMP::RunImpl() {
  if (input_data_.empty()) {
    result_data_.clear();
    return true;
  }

  auto blocks = MakeSortedBlocks(input_data_);
  result_data_ = MergeBlocks(std::move(blocks));
  return true;
}

bool DoubleSortEvenOddBatcherOMP::PostProcessingImpl() {
  GetOutput() = result_data_;
  return true;
}

}  // namespace gusev_d_double_sort_even_odd_batcher_omp_task_threads
