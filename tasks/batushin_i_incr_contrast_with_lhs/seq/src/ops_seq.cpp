#include "batushin_i_incr_contrast_with_lhs/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "batushin_i_incr_contrast_with_lhs/common/include/common.hpp"

namespace batushin_i_incr_contrast_with_lhs {

BatushinITestTaskSEQ::BatushinITestTaskSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().resize(in.size());
}

bool BatushinITestTaskSEQ::ValidationImpl() {
  return !GetInput().empty();
}

bool BatushinITestTaskSEQ::PreProcessingImpl() {
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

std::pair<unsigned char, unsigned char> FindMinMax(const std::vector<unsigned char> &data) {
  if (data.empty()) {
    return {0, 0};
  }

  unsigned char minimum = data[0];
  unsigned char maximum = data[0];

  for (size_t idx = 1; idx < data.size(); ++idx) {
    minimum = std::min(minimum, data[idx]);
    maximum = std::max(maximum, data[idx]);
  }

  return {minimum, maximum};
}

void FillUniformImage(std::vector<unsigned char> &output, size_t size) {
  output.assign(size, 128);
}

}  // namespace

bool BatushinITestTaskSEQ::RunImpl() {
  const std::vector<unsigned char> &source = GetInput();
  std::vector<unsigned char> &destination = GetOutput();

  auto [min_value, max_value] = FindMinMax(source);

  if (min_value == max_value) {
    FillUniformImage(destination, source.size());
    return true;
  }

  const double scale_coefficient = 255.0 / static_cast<double>(max_value - min_value);

  destination.resize(source.size());
  for (size_t position = 0; position < source.size(); ++position) {
    destination[position] = NormalizePixel(source[position], min_value, scale_coefficient);
  }

  return true;
}

bool BatushinITestTaskSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace batushin_i_incr_contrast_with_lhs
