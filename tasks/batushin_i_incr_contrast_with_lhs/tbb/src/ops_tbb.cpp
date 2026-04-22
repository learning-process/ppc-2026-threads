#include "batushin_i_incr_contrast_with_lhs/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "batushin_i_incr_contrast_with_lhs/common/include/common.hpp"
#include "oneapi/tbb/parallel_for.h"

namespace batushin_i_incr_contrast_with_lhs {

BatushinIIncrContrastWithLhsTBB::BatushinIIncrContrastWithLhsTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().resize(in.size());
}

bool BatushinIIncrContrastWithLhsTBB::ValidationImpl() {
  return !GetInput().empty();
}

bool BatushinIIncrContrastWithLhsTBB::PreProcessingImpl() {
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

  auto result = tbb::parallel_reduce(
      tbb::blocked_range<std::size_t>(0, data.size()), std::make_pair(data[0], data[0]),
      [&](const tbb::blocked_range<std::size_t> &r, std::pair<unsigned char, unsigned char> local) {
    for (std::size_t i = r.begin(); i != r.end(); ++i) {
      local.first = std::min(local.first, data[i]);
      local.second = std::max(local.second, data[i]);
    }
    return local;
  }, [](const std::pair<unsigned char, unsigned char> &a, const std::pair<unsigned char, unsigned char> &b) {
    return std::make_pair(std::min(a.first, b.first), std::max(a.second, b.second));
  });

  return result;
}

void NormalizeImage(const std::vector<unsigned char> &source, std::vector<unsigned char> &destination,
                    unsigned char min_value, double scale_coefficient) {
  tbb::parallel_for(tbb::blocked_range<std::size_t>(0, source.size()),
                    [&](const tbb::blocked_range<std::size_t> &range) {
    for (std::size_t i = range.begin(); i != range.end(); ++i) {
      destination[i] = NormalizePixel(source[i], min_value, scale_coefficient);
    }
  });
}

void FillUniformImage(std::vector<unsigned char> &output, std::size_t size) {
  tbb::parallel_for(tbb::blocked_range<std::size_t>(0, size), [&](const tbb::blocked_range<std::size_t> &range) {
    for (std::size_t i = range.begin(); i != range.end(); ++i) {
      output[i] = 128;
    }
  });
}

}  // namespace

bool BatushinIIncrContrastWithLhsTBB::RunImpl() {
  const std::vector<unsigned char> &source = GetInput();
  std::vector<unsigned char> &destination = GetOutput();

  auto min_max = FindMinMaxParallel(source);
  unsigned char min_value = min_max.first;
  unsigned char max_value = min_max.second;

  if (min_value == max_value) {
    FillUniformImage(destination, source.size());
    return true;
  }

  const double scale_coefficient = 255.0 / static_cast<double>(max_value - min_value);
  destination.resize(source.size());

  NormalizeImage(source, destination, min_value, scale_coefficient);

  return true;
}

bool BatushinIIncrContrastWithLhsTBB::PostProcessingImpl() {
  return true;
}

}  // namespace batushin_i_incr_contrast_with_lhs
