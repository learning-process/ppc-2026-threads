#include "klimenko_v_lsh_contrast_incr_omp/omp/include/ops_omp.hpp"

#include <atomic>
#include <numeric>
#include <vector>

#include "klimenko_v_lsh_contrast_incr_omp/common/include/common.hpp"
#include "util/include/util.hpp"

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

  const int size = static_cast<int>(input.size());

  int min_val = input[0];
  int max_val = input[0];

#pragma omp parallel for reduction(min : min_val) reduction(max : max_val)
  for (int i = 0; i < size; ++i) {
    if (input[i] < min_val) {
      min_val = input[i];
    }
    if (input[i] > max_val) {
      max_val = input[i];
    }
  }

  if (max_val == min_val) {
#pragma omp parallel for
    for (int i = 0; i < size; ++i) {
      output[i] = input[i];
    }
    return true;
  }

#pragma omp parallel for
  for (int i = 0; i < static_cast<int>(input.size()); ++i) {
    output[i] = ((input[i] - min_val) * 255) / (max_val - min_val);
  }
  return true;
}

bool KlimenkoVLSHContrastIncrOMP::PostProcessingImpl() {
  return GetOutput().size() == GetInput().size();
}

}  // namespace klimenko_v_lsh_contrast_incr_omp
