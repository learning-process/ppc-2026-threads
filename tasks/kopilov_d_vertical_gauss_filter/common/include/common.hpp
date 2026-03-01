#pragma once

#include <cstdint>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace kopilov_d_vertical_gauss_filter {

struct Matrix {
  int width = 0;
  int height = 0;
  std::vector<std::uint8_t> data;
};

using InType = Matrix;
using OutType = Matrix;
using TestType = std::tuple<InType, OutType>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace kopilov_d_vertical_gauss_filter