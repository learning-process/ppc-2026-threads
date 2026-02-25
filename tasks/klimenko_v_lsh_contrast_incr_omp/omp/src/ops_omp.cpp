#include "klimenko_v_lsh_contrast_incr_omp/omp/include/ops_omp.hpp"

#include <algorithm>
#include <cstddef>
#include <vector>

#include "klimenko_v_lsh_contrast_incr_omp/common/include/common.hpp"

namespace klimenko_v_lsh_contrast_incr_omp {

KlimenkoVLSHContrastIncrOMP::KlimenkoVLSHContrastIncrOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool KlimenkoVLSHContrastIncrOMP::ValidationImpl() {
  return !GetInput().empty();
}

bool KlimenkoVLSHContrastIncrOMP::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return true;
}

bool KlimenkoVLSHContrastIncrOMP::RunImpl() {
  const auto &input = GetInput();
  auto &output = GetOutput();

  if (input.empty()) {
    return false;
  }

  const size_t size = input.size();

  int min_val = input[0];
  int max_val = input[0];

#pragma omp parallel for default(none) shared(input, size) reduction(min : min_val) reduction(max : max_val)
  for (size_t i = 0; i < size; ++i) {
    min_val = std::min(input[i], min_val);
    max_val = std::max(input[i], max_val);
  }

  if (max_val == min_val) {
#pragma omp parallel for default(none) shared(input, output, size)
    for (size_t i = 0; i < size; ++i) {
      output[i] = input[i];
    }
    return true;
  }
  const int range = max_val - min_val;
#pragma omp parallel for default(none) shared(input, output, size, min_val, range)
  for (size_t i = 0; i < size; ++i) {
    output[i] = ((input[i] - min_val) * 255) / range;
  }
  return true;
}

bool KlimenkoVLSHContrastIncrOMP::PostProcessingImpl() {
  return GetOutput().size() == GetInput().size();
}

}  // namespace klimenko_v_lsh_contrast_incr_omp
