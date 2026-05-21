#include "batushin_i_incr_contrast_with_lhs/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <thread>
#include <utility>
#include <vector>

#include "batushin_i_incr_contrast_with_lhs/common/include/common.hpp"

namespace batushin_i_incr_contrast_with_lhs {

BatushinIIncrContrastWithLhsSTL::BatushinIIncrContrastWithLhsSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().resize(in.size());
}

bool BatushinIIncrContrastWithLhsSTL::ValidationImpl() {
  return !GetInput().empty();
}

bool BatushinIIncrContrastWithLhsSTL::PreProcessingImpl() {
  return true;
}

namespace {

unsigned char NormalizePixel(unsigned char pixel, unsigned char min_val, double scale_factor) {
  double normalized = static_cast<double>(pixel - min_val) * scale_factor;
  normalized = std::floor(normalized + 0.5);
  normalized = std::max(normalized, 0.0);
  normalized = std::min(normalized, 255.0);
  return static_cast<unsigned char>(normalized);
}

unsigned int GetNumThreads() {
  unsigned int threads = std::thread::hardware_concurrency();
  return (threads == 0) ? 1 : threads;
}

size_t ComputeChunkSize(size_t data_size, unsigned int num_threads) {
  return std::max(static_cast<size_t>(1), data_size / num_threads);
}

size_t ComputeNumBlocks(size_t data_size, size_t chunk_size) {
  return (data_size + chunk_size - 1) / chunk_size;
}

struct MinMaxResult {
  unsigned char min_val;
  unsigned char max_val;
};

MinMaxResult FindMinMaxInBlock(const std::vector<unsigned char> &data, size_t start, size_t end) {
  unsigned char local_min = 255;
  unsigned char local_max = 0;
  for (size_t idx = start; idx < end; ++idx) {
    unsigned char val = data[idx];
    local_min = std::min(local_min, val);
    local_max = std::max(local_max, val);
  }
  // Используем designated initializers
  return MinMaxResult{.min_val = local_min, .max_val = local_max};
}

MinMaxResult MergeMinMaxResults(const std::vector<MinMaxResult> &results) {
  unsigned char global_min = 255;
  unsigned char global_max = 0;
  for (const auto &result : results) {
    global_min = std::min(global_min, result.min_val);
    global_max = std::max(global_max, result.max_val);
  }
  // Используем designated initializers
  return MinMaxResult{.min_val = global_min, .max_val = global_max};
}

std::pair<unsigned char, unsigned char> FindMinMaxParallel(const std::vector<unsigned char> &data) {
  if (data.empty()) {
    return {0, 0};
  }

  const size_t data_size = data.size();
  const unsigned int num_threads = GetNumThreads();
  const size_t chunk_size = ComputeChunkSize(data_size, num_threads);
  const size_t num_blocks = ComputeNumBlocks(data_size, chunk_size);

  std::vector<MinMaxResult> block_results(num_blocks, MinMaxResult{.min_val = 255, .max_val = 0});
  std::vector<std::thread> threads;
  threads.reserve(num_blocks);

  for (size_t block_idx = 0; block_idx < num_blocks; ++block_idx) {
    size_t start = block_idx * chunk_size;
    size_t end = std::min(start + chunk_size, data_size);

    threads.emplace_back([&data, &block_results, block_idx, start, end] {
      block_results[block_idx] = FindMinMaxInBlock(data, start, end);
    });
  }

  for (auto &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }

  MinMaxResult final_result = MergeMinMaxResults(block_results);
  return {final_result.min_val, final_result.max_val};
}

void NormalizeBlock(const std::vector<unsigned char> &source, std::vector<unsigned char> &destination, size_t start,
                    size_t end, unsigned char min_value, double scale_coefficient) {
  for (size_t idx = start; idx < end; ++idx) {
    destination[idx] = NormalizePixel(source[idx], min_value, scale_coefficient);
  }
}

void NormalizeImageParallel(const std::vector<unsigned char> &source, std::vector<unsigned char> &destination,
                            unsigned char min_value, double scale_coefficient) {
  const size_t data_size = source.size();
  const unsigned int num_threads = GetNumThreads();
  const size_t chunk_size = ComputeChunkSize(data_size, num_threads);
  const size_t num_blocks = ComputeNumBlocks(data_size, chunk_size);

  std::vector<std::thread> threads;
  threads.reserve(num_blocks);

  for (size_t block_idx = 0; block_idx < num_blocks; ++block_idx) {
    size_t start = block_idx * chunk_size;
    size_t end = std::min(start + chunk_size, data_size);

    threads.emplace_back([&source, &destination, start, end, min_value, scale_coefficient] {
      NormalizeBlock(source, destination, start, end, min_value, scale_coefficient);
    });
  }

  for (auto &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }
}

void FillUniformBlock(std::vector<unsigned char> &output, size_t start, size_t end) {
  for (size_t idx = start; idx < end; ++idx) {
    output[idx] = 128;
  }
}

void FillUniformImageParallel(std::vector<unsigned char> &output, size_t size) {
  const unsigned int num_threads = GetNumThreads();
  const size_t chunk_size = ComputeChunkSize(size, num_threads);
  const size_t num_blocks = ComputeNumBlocks(size, chunk_size);

  std::vector<std::thread> threads;
  threads.reserve(num_blocks);

  for (size_t block_idx = 0; block_idx < num_blocks; ++block_idx) {
    size_t start = block_idx * chunk_size;
    size_t end = std::min(start + chunk_size, size);

    threads.emplace_back([&output, start, end] { FillUniformBlock(output, start, end); });
  }

  for (auto &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }
}

}  // namespace

bool BatushinIIncrContrastWithLhsSTL::RunImpl() {
  const std::vector<unsigned char> &source = GetInput();
  std::vector<unsigned char> &destination = GetOutput();

  auto min_max = FindMinMaxParallel(source);
  unsigned char min_value = min_max.first;
  unsigned char max_value = min_max.second;

  if (min_value == max_value) {
    FillUniformImageParallel(destination, source.size());
    return true;
  }

  const double scale_coefficient = 255.0 / static_cast<double>(max_value - min_value);
  destination.resize(source.size());

  NormalizeImageParallel(source, destination, min_value, scale_coefficient);

  return true;
}

bool BatushinIIncrContrastWithLhsSTL::PostProcessingImpl() {
  return true;
}

}  // namespace batushin_i_incr_contrast_with_lhs
