#include "lukin_i_ench_contr_lin_hist/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <algorithm>
#include <utility>
#include <vector>

#include "lukin_i_ench_contr_lin_hist/common/include/common.hpp"
#include "oneapi/tbb/parallel_for.h"

namespace lukin_i_ench_contr_lin_hist {

using MinMax = std::pair<unsigned char, unsigned char>;

LukinITestTaskTBB::LukinITestTaskTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType(GetInput().size());
}

bool LukinITestTaskTBB::ValidationImpl() {
  return !(GetInput().empty());
}

bool LukinITestTaskTBB::PreProcessingImpl() {
  return true;
}

bool LukinITestTaskTBB::RunImpl() {
  const auto &input = GetInput();
  auto &output = GetOutput();

  unsigned char min = 255;
  unsigned char max = 0;

  const int size = static_cast<int>(input.size());

  auto result = tbb::parallel_reduce(tbb::blocked_range<int>(0, size), MinMax(min, max),
                                     [&](const tbb::blocked_range<int> &r, MinMax local) {
    unsigned char local_min = local.first;
    unsigned char local_max = local.second;

    for (int i = r.begin(); i < r.end(); i++) {
      unsigned char curr = input[i];
      local_max = std::max(local_max, curr);
      local_min = std::min(local_min, curr);
    }

    return MinMax(local_min, local_max);
  }, [](MinMax left, MinMax right) {
    return MinMax(std::min(left.first, right.first), std::max(left.second, right.second));
  });

  min = result.first;
  max = result.second;

  if (max == min) {
    output = input;
    return true;
  }

  float scale = 255.0F / static_cast<float>(max - min);

  tbb::parallel_for(tbb::blocked_range<int>(0, size), [&](const tbb::blocked_range<int> &r) {
    for (int i = r.begin(); i < r.end(); i++) {
      output[i] = static_cast<unsigned char>(static_cast<float>(input[i] - min) * scale);
    }
  });

  return true;
}

bool LukinITestTaskTBB::PostProcessingImpl() {
  return true;
}

}  // namespace lukin_i_ench_contr_lin_hist
