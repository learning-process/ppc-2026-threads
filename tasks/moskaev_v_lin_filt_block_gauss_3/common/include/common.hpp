#pragma once

#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace moskaev_v_lin_filt_block_gauss_3 {

using InType = std::tuple<int, int, int, int, std::vector<uint8_t>>;
using OutType = std::vector<uint8_t>;
using TestType = std::tuple<int, int, int, int, std::vector<uint8_t>, std::vector<uint8_t>, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

const std::vector<float> kGaussianKernel = {1.0F / 16, 2.0F / 16, 1.0F / 16, 2.0F / 16, 4.0F / 16,
                                            2.0F / 16, 1.0F / 16, 2.0F / 16, 1.0F / 16};

struct ImageInfo {
  int width{0};
  int height{0};
  int channels{0};
  std::vector<uint8_t> data;
};

}  // namespace moskaev_v_lin_filt_block_gauss_3
