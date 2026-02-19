#pragma once

#include <string>
#include <tuple>

#include "crs.hpp"
#include "task/include/task.hpp"

namespace guseva_crs {

struct CRS {
  size_t nz{};
  size_t nrows{};
  size_t ncols{};
  std::vector<double> values;
  std::vector<size_t> cols;
  std::vector<size_t> row_ptrs;
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
