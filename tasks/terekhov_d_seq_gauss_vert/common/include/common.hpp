#pragma once

#include <vector>

#include "task/include/task.hpp"

namespace terekhov_d_seq_gauss_vert {

struct Image {
  int width = 0;
  int height = 0;
  std::vector<int> data;
};

using InType = Image;
using OutType = Image;
using TestType = int;
using BaseTask = ppc::task::Task<InType, OutType>;

const std::vector<float> kGaussKernel = {1.0F / 16, 2.0F / 16, 1.0F / 16, 2.0F / 16, 4.0F / 16,
                                         2.0F / 16, 1.0F / 16, 2.0F / 16, 1.0F / 16};

}  // namespace terekhov_d_seq_gauss_vert
