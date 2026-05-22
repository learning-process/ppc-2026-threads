#include "klimenko_v_lsh_contrast_incr/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <algorithm>
#include <cstddef>
#include <vector>

#include "klimenko_v_lsh_contrast_incr/common/include/common.hpp"
#include "oneapi/tbb/parallel_for.h"

namespace klimenko_v_lsh_contrast_incr {

KlimenkoVLSHContrastIncrTBB::KlimenkoVLSHContrastIncrTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool KlimenkoVLSHContrastIncrTBB::ValidationImpl() {
  return !GetInput().empty();
}

bool KlimenkoVLSHContrastIncrTBB::PreProcessingImpl() {
  GetOutput().resize(GetInput().size());
  return true;
}

bool KlimenkoVLSHContrastIncrTBB::RunImpl() {
  const auto &input = GetInput();
  auto &output = GetOutput();

  if (input.empty()) {
    return false;
  }

  const size_t size = input.size();

  struct MinMax {
    int min_val;
    int max_val;
  };

  MinMax result =
      tbb::parallel_reduce(tbb::blocked_range<size_t>(0, size), MinMax{.min_val = input[0], .max_val = input[0]},
                           [&](const tbb::blocked_range<size_t> &r, MinMax local) {
    for (size_t i = r.begin(); i < r.end(); ++i) {
      local.min_val = std::min(local.min_val, input[i]);
      local.max_val = std::max(local.max_val, input[i]);
    }
    return local;
  }, [](const MinMax &a, const MinMax &b) {
    return MinMax{.min_val = std::min(a.min_val, b.min_val), .max_val = std::max(a.max_val, b.max_val)};
  });

  int min_val = result.min_val;
  int max_val = result.max_val;

  if (max_val == min_val) {
    tbb::parallel_for(tbb::blocked_range<size_t>(0, size), [&](const tbb::blocked_range<size_t> &r) {
      for (size_t i = r.begin(); i < r.end(); ++i) {
        output[i] = input[i];
      }
    });
    return true;
  }

  tbb::parallel_for(tbb::blocked_range<size_t>(0, size), [&](const tbb::blocked_range<size_t> &r) {
    for (size_t i = r.begin(); i < r.end(); ++i) {
      output[i] = ((input[i] - min_val) * 255) / (max_val - min_val);
    }
  });

  return true;
}

bool KlimenkoVLSHContrastIncrTBB::PostProcessingImpl() {
  return GetOutput().size() == GetInput().size();
}

}  // namespace klimenko_v_lsh_contrast_incr
