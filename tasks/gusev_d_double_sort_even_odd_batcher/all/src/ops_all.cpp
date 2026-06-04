#include "gusev_d_double_sort_even_odd_batcher/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

#include "gusev_d_double_sort_even_odd_batcher/all/include/common.hpp"
#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/global_control.h"
#include "oneapi/tbb/parallel_for.h"
#include "util/include/util.hpp"

namespace gusev_d_double_sort_even_odd_batcher_all_task_threads {
namespace {

constexpr int kRadixPasses = 8;
constexpr int kBitsPerByte = 8;
constexpr size_t kRadixBuckets = 256;
constexpr uint64_t kBucketMask = 0xFFULL;
constexpr size_t kMinParallelElements = 128;
constexpr size_t kMinThreadedTasks = 2;
constexpr size_t kMinTbbComparePairs = 512;

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

struct MpiContext {
  bool active = false;
  int rank = 0;
  int size = 1;
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
      count.at(GetBucketIndex(value, shift))++;
    }
    BuildPrefixSums(count);

    for (ValueType value : *source) {
      const auto bucket = GetBucketIndex(value, shift);
      (*destination)[count.at(bucket)++] = value;
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

void CompareExchangeBlocks(std::vector<MergeItem> &data, size_t first, size_t distance) {
  if (distance < kMinTbbComparePairs) {
    for (size_t i = 0; i < distance; ++i) {
      CompareExchange(data, first + i, first + distance + i);
    }
    return;
  }

  tbb::parallel_for(tbb::blocked_range<size_t>(0, distance), [&](const tbb::blocked_range<size_t> &range) {
    for (size_t i = range.begin(); i != range.end(); ++i) {
      CompareExchange(data, first + i, first + distance + i);
    }
  });
}

void OddEvenMerge(std::vector<MergeItem> &data) {
  if (data.size() < 2) {
    return;
  }

  auto distance = data.size() / 2;
  CompareExchangeBlocks(data, 0, distance);

  for (distance /= 2; distance > 0; distance /= 2) {
    const auto step = distance * 2;
    for (size_t first = distance; (first + distance) < data.size(); first += step) {
      CompareExchangeBlocks(data, first, distance);
    }
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
  OddEvenMerge(items);
  return ExtractMergedValues(items, left.size() + right.size());
}

BlockRange GetBlockRange(size_t block_index, size_t block_count, size_t total_size) {
  return {
      .begin = (block_index * total_size) / block_count,
      .end = ((block_index + 1) * total_size) / block_count,
  };
}

size_t GetBlockCount(size_t input_size, size_t parallelism) {
  const auto target_blocks = std::max<size_t>(1, parallelism * 2);
  return std::max<size_t>(1, std::min(input_size, target_blocks));
}

void FillAndSortBlock(const InType &input, Block &block, BlockRange range) {
  block.assign(input.begin() + static_cast<std::ptrdiff_t>(range.begin),
               input.begin() + static_cast<std::ptrdiff_t>(range.end));
  RadixSortDoubles(block);
}

void StoreCurrentException(std::exception_ptr &worker_exception, std::mutex &exception_mutex) noexcept {
  const std::scoped_lock lock(exception_mutex);
  if (worker_exception == nullptr) {
    worker_exception = std::current_exception();
  }
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

  const auto signed_block_count = static_cast<std::ptrdiff_t>(block_count);
  const auto omp_threads = static_cast<int>(std::min<size_t>(parallelism, block_count));
  static_cast<void>(omp_threads);

  std::exception_ptr worker_exception;
  std::mutex exception_mutex;
#pragma omp parallel for schedule(static) num_threads(omp_threads) default(none) \
    shared(input, blocks, block_count, total_size, signed_block_count, worker_exception, exception_mutex)
  for (std::ptrdiff_t block = 0; block < signed_block_count; ++block) {
    try {
      const auto block_index = static_cast<size_t>(block);
      FillAndSortBlock(input, blocks[block_index], GetBlockRange(block_index, block_count, total_size));
    } catch (...) {
      StoreCurrentException(worker_exception, exception_mutex);
    }
  }

  if (worker_exception != nullptr) {
    std::rethrow_exception(worker_exception);
  }

  return blocks;
}

template <typename Function>
void RunSequentialByIndex(size_t work_count, const Function &function) {
  for (size_t index = 0; index < work_count; ++index) {
    function(index);
  }
}

template <typename Function>
void RunThreadChunk(size_t thread_index, size_t work_count, size_t thread_count, const Function &function,
                    std::exception_ptr &worker_exception, std::mutex &exception_mutex) noexcept {
  try {
    for (size_t index = thread_index; index < work_count; index += thread_count) {
      function(index);
    }
  } catch (...) {
    StoreCurrentException(worker_exception, exception_mutex);
  }
}

void JoinWorkersAndRethrow(std::vector<std::thread> &workers, const std::exception_ptr &worker_exception) {
  for (auto &worker : workers) {
    worker.join();
  }

  if (worker_exception != nullptr) {
    std::rethrow_exception(worker_exception);
  }
}

template <typename Function>
void RunThreadedByIndex(size_t work_count, size_t max_threads, const Function &function) {
  if (work_count < kMinThreadedTasks || max_threads < kMinThreadedTasks) {
    RunSequentialByIndex(work_count, function);
    return;
  }

  const auto thread_count = std::min(work_count, max_threads);
  std::vector<std::thread> workers;
  workers.reserve(thread_count);
  std::exception_ptr worker_exception;
  std::mutex exception_mutex;
  for (size_t thread_index = 0; thread_index < thread_count; ++thread_index) {
    workers.emplace_back([&, thread_index] {
      RunThreadChunk(thread_index, work_count, thread_count, function, worker_exception, exception_mutex);
    });
  }

  JoinWorkersAndRethrow(workers, worker_exception);
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

Block SortLocal(const InType &input, size_t parallelism) {
  if (input.empty()) {
    return {};
  }

  auto blocks = MakeSortedBlocks(input, parallelism);
  return MergeBlocks(std::move(blocks), parallelism);
}

int ToMpiCount(size_t size) {
  if (size > static_cast<size_t>(std::numeric_limits<int>::max())) {
    throw std::runtime_error("Input is too large for MPI int counts");
  }
  return static_cast<int>(size);
}

MpiContext GetMpiContext() {
  int initialized = 0;
  MPI_Initialized(&initialized);
  if (initialized == 0) {
    return {};
  }

  int finalized = 0;
  MPI_Finalized(&finalized);
  if (finalized != 0) {
    return {};
  }

  MpiContext context{.active = true};
  MPI_Comm_rank(MPI_COMM_WORLD, &context.rank);
  MPI_Comm_size(MPI_COMM_WORLD, &context.size);
  return context;
}

std::vector<int> BuildCounts(size_t total_size, int parts) {
  std::vector<int> counts(static_cast<size_t>(parts));
  const auto base = total_size / static_cast<size_t>(parts);
  const auto remainder = total_size % static_cast<size_t>(parts);

  for (int part = 0; part < parts; ++part) {
    const auto part_size = base + (std::cmp_less(part, remainder) ? 1 : 0);
    counts[static_cast<size_t>(part)] = ToMpiCount(part_size);
  }

  return counts;
}

std::vector<int> BuildDisplacements(const std::vector<int> &counts) {
  std::vector<int> displacements(counts.size());
  for (size_t i = 1; i < counts.size(); ++i) {
    displacements[i] = displacements[i - 1] + counts[i - 1];
  }
  return displacements;
}

Block ScatterInput(const InType &input, MpiContext context) {
  if (!context.active || context.size == 1) {
    return input;
  }

  uint64_t total_size_wire = context.rank == 0 ? static_cast<uint64_t>(input.size()) : 0;
  MPI_Bcast(&total_size_wire, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
  const auto total_size = static_cast<size_t>(total_size_wire);

  const auto counts = BuildCounts(total_size, context.size);
  const auto displacements = BuildDisplacements(counts);
  Block local_data(static_cast<size_t>(counts[static_cast<size_t>(context.rank)]));

  MPI_Scatterv(context.rank == 0 ? input.data() : nullptr, counts.data(), displacements.data(), MPI_DOUBLE,
               local_data.data(), counts[static_cast<size_t>(context.rank)], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  return local_data;
}

BlockList SplitGatheredBlocks(const Block &gathered_data, const std::vector<int> &counts) {
  BlockList blocks;
  blocks.reserve(counts.size());

  size_t offset = 0;
  for (int count : counts) {
    const auto block_size = static_cast<size_t>(count);
    blocks.emplace_back(gathered_data.begin() + static_cast<std::ptrdiff_t>(offset),
                        gathered_data.begin() + static_cast<std::ptrdiff_t>(offset + block_size));
    offset += block_size;
  }

  return blocks;
}

Block BroadcastResult(Block result, MpiContext context) {
  int result_size = context.rank == 0 ? ToMpiCount(result.size()) : 0;
  MPI_Bcast(&result_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
  result.resize(static_cast<size_t>(result_size));
  MPI_Bcast(result.data(), result_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  return result;
}

Block GatherAndMerge(const Block &local_sorted, MpiContext context, size_t parallelism) {
  if (!context.active || context.size == 1) {
    return local_sorted;
  }

  const auto local_count = ToMpiCount(local_sorted.size());
  std::vector<int> counts(static_cast<size_t>(context.size));
  MPI_Gather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

  Block gathered_data;
  std::vector<int> displacements;
  if (context.rank == 0) {
    displacements = BuildDisplacements(counts);
    const auto total_count = displacements.back() + counts.back();
    gathered_data.resize(static_cast<size_t>(total_count));
  }

  MPI_Gatherv(local_sorted.data(), local_count, MPI_DOUBLE, context.rank == 0 ? gathered_data.data() : nullptr,
              counts.data(), context.rank == 0 ? displacements.data() : nullptr, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  Block result;
  if (context.rank == 0) {
    result = MergeBlocks(SplitGatheredBlocks(gathered_data, counts), parallelism);
  }

  return BroadcastResult(std::move(result), context);
}

}  // namespace

DoubleSortEvenOddBatcherALL::DoubleSortEvenOddBatcherALL(const InType &in) : input_data_(in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

bool DoubleSortEvenOddBatcherALL::ValidationImpl() {
  return GetOutput().empty();
}

bool DoubleSortEvenOddBatcherALL::PreProcessingImpl() {
  input_data_ = GetInput();
  result_data_.clear();
  return true;
}

bool DoubleSortEvenOddBatcherALL::RunImpl() {
  const auto parallelism = static_cast<size_t>(std::max(1, ppc::util::GetNumThreads()));
  const tbb::global_control control(tbb::global_control::max_allowed_parallelism, parallelism);
  const auto mpi_context = GetMpiContext();

  const auto local_data = ScatterInput(input_data_, mpi_context);
  const auto local_sorted = SortLocal(local_data, parallelism);
  result_data_ = GatherAndMerge(local_sorted, mpi_context, parallelism);
  return true;
}

bool DoubleSortEvenOddBatcherALL::PostProcessingImpl() {
  GetOutput() = result_data_;
  return true;
}

}  // namespace gusev_d_double_sort_even_odd_batcher_all_task_threads
