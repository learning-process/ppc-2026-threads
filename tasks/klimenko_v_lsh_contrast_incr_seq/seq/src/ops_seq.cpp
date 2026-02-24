#include "klimenko_v_lsh_contrast_incr_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <vector>

#include "klimenko_v_lsh_contrast_incr_seq/common/include/common.hpp"

namespace klimenko_v_lsh_contrast_incr_seq {

KlimenkoVLSHContrastIncrSEQ::KlimenkoVLSHContrastIncrSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool KlimenkoVLSHContrastIncrSEQ::ValidationImpl() {
  return !GetInput().empty();
}

bool KlimenkoVLSHContrastIncrSEQ::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return true;
}

bool KlimenkoVLSHContrastIncrSEQ::RunImpl() {
  const auto &input = GetInput();
  auto &output = GetOutput();

  if (input.empty()) {
    return false;
  }

  auto minmax = std::ranges::minmax_element(input);
  int min_val = *minmax.min;
  int max_val = *minmax.max;

  if (max_val == min_val) {
    output = input;
    return true;
  }

  for (size_t i = 0; i < input.size(); ++i) {
    output[i] = ((input[i] - min_val) * 255) / (max_val - min_val);
  }

  return true;
}

bool KlimenkoVLSHContrastIncrSEQ::PostProcessingImpl() {
  return GetOutput().size() == GetInput().size();
}

}  // namespace klimenko_v_lsh_contrast_incr_seq
