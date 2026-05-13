#include "batushin_i_incr_contrast_with_lhs/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <mutex>
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

std::pair<unsigned char, unsigned char> FindMinMaxParallel(const std::vector<unsigned char> &data) {
  if (data.empty()) {
    return {0, 0};
  }

  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 1;
  }

  const size_t data_size = data.size();
  const size_t chunk_size = std::max(static_cast<size_t>(1), data_size / num_threads);
  const size_t num_blocks = (data_size + chunk_size - 1) / chunk_size;

  std::vector<std::pair<unsigned char, unsigned char>> block_results(num_blocks, {255, 0});
  std::vector<std::thread> threads;
  threads.reserve(num_blocks);

  for (size_t b = 0; b < num_blocks; ++b) {
    size_t start = b * chunk_size;
    size_t end = std::min(start + chunk_size, data_size);

    threads.emplace_back([&data, &block_results, b, start, end] {
      unsigned char local_min = 255;
      unsigned char local_max = 0;
      for (size_t i = start; i < end; ++i) {
        unsigned char val = data[i];
        if (val < local_min) {
          local_min = val;
        }
        if (val > local_max) {
          local_max = val;
        }
      }
      block_results[b] = {local_min, local_max};
    });
  }

  for (auto &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }

  unsigned char global_min = 255;
  unsigned char global_max = 0;
  for (const auto &res : block_results) {
    if (res.first < global_min) {
      global_min = res.first;
    }
    if (res.second > global_max) {
      global_max = res.second;
    }
  }

  return {global_min, global_max};
}

void NormalizeImageParallel(const std::vector<unsigned char> &source, std::vector<unsigned char> &destination,
                            unsigned char min_value, double scale_coefficient) {
  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 1;
  }

  const size_t data_size = source.size();
  const size_t chunk_size = std::max(static_cast<size_t>(1), data_size / num_threads);
  const size_t num_blocks = (data_size + chunk_size - 1) / chunk_size;

  std::vector<std::thread> threads;
  threads.reserve(num_blocks);

  for (size_t b = 0; b < num_blocks; ++b) {
    size_t start = b * chunk_size;
    size_t end = std::min(start + chunk_size, data_size);

    threads.emplace_back([&source, &destination, start, end, min_value, scale_coefficient] {
      for (size_t i = start; i < end; ++i) {
        destination[i] = NormalizePixel(source[i], min_value, scale_coefficient);
      }
    });
  }

  for (auto &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }
}

void FillUniformImageParallel(std::vector<unsigned char> &output, size_t size) {
  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 1;
  }

  const size_t chunk_size = std::max(static_cast<size_t>(1), size / num_threads);
  const size_t num_blocks = (size + chunk_size - 1) / chunk_size;

  std::vector<std::thread> threads;
  threads.reserve(num_blocks);

  for (size_t b = 0; b < num_blocks; ++b) {
    size_t start = b * chunk_size;
    size_t end = std::min(start + chunk_size, size);

    threads.emplace_back([&output, start, end] {
      for (size_t i = start; i < end; ++i) {
        output[i] = 128;
      }
    });
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
