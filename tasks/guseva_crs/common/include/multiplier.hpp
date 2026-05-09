#pragma once

#include <cstddef>

#include "common.hpp"
namespace guseva_crs {

class Multiplier {
 public:
  Multiplier() = default;
  virtual ~Multiplier() = default;

  [[nodiscard]] virtual CRS Transpose(const CRS &a) const {
    CRS transposed;
    transposed.nz = a.nz;
    transposed.ncols = a.nrows;
    transposed.nrows = a.ncols;
    transposed.row_ptrs.resize(transposed.nrows + 1, 0);
    transposed.values.resize(a.nz, 0);
    transposed.cols.resize(a.nz, 0);
    for (std::size_t i = 0; i < a.nz; i++) {
      transposed.row_ptrs[a.cols[i] + 1]++;
    }
    std::size_t row = 0;
    for (std::size_t i = 1; i <= transposed.nrows; i++) {
      auto tmp = transposed.row_ptrs[i];
      transposed.row_ptrs[i] = row;
      row += tmp;
    }
    for (std::size_t i = 0; i < a.nrows; i++) {
      auto j1 = a.row_ptrs[i];
      auto j2 = a.row_ptrs[i + 1];
      auto col = i;
      for (auto j = j1; j < j2; j++) {
        auto value = a.values[j];
        auto row_ind = a.cols[j];
        auto ind = transposed.row_ptrs[row_ind + 1];
        transposed.values[ind] = value;
        transposed.cols[ind] = col;
        transposed.row_ptrs[row_ind + 1]++;
      }
    }
    return transposed;
  };

  [[nodiscard]] virtual CRS Multiply(const CRS &a, const CRS &b) const = 0;
};

}  // namespace guseva_crs
