#pragma once

#include <vector>

#include "task/include/task.hpp"

namespace terekhov_d_seq_gauss_vert {

struct Image {
  int width;
  int height;
  std::vector<int> data;
};

using InType = Image;
using OutType = Image;
using TestType = int;
using BaseTask = ppc::task::Task<InType, OutType>;

// Ядро Гаусса 3x3
const std::vector<float> kGaussKernel = {
    1.0f / 16, 2.0f / 16, 1.0f / 16,
    2.0f / 16, 4.0f / 16, 2.0f / 16,
    1.0f / 16, 2.0f / 16, 1.0f / 16
};

}  // namespace terekhov_d_seq_gauss_vert
