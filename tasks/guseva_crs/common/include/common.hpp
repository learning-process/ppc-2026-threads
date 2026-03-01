#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace guseva_crs {

struct CRS {
  std::size_t nz{};
  std::size_t nrows{};
  std::size_t ncols{};
  std::vector<double> values;
  std::vector<std::size_t> cols;
  std::vector<std::size_t> row_ptrs;
};

using InType = std::tuple<CRS, CRS>;
using OutType = CRS;
using TestType = std::string;
using BaseTask = ppc::task::Task<InType, OutType>;

constexpr double kZERO = 10e-5;

inline bool Equal(const CRS &a, const CRS &b) {
  return a.nz == b.nz && a.ncols == b.ncols && a.nrows == b.nrows &&
         std::ranges::equal(a.values, b.values, [](double a, double b) { return std::fabs(a - b) < kZERO; });
}

}  // namespace guseva_crs
