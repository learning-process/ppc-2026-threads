#pragma once
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "guseva_crs/common/include/common.hpp"
#include "guseva_crs/common/include/multiplier.hpp"

namespace guseva_crs {

class MultiplierSeq : public Multiplier {
 public:
  [[nodiscard]] CRS Multiply(const CRS &a, const CRS &b) const override {
    CRS result;
    result.nrows = a.nrows;
    result.ncols = b.ncols;
    result.row_ptrs.push_back(0);

    auto bt = this->Transpose(b);
    double sum = 0;
    std::size_t nz = 0;

    std::vector<int> temp(a.nrows, -1);
    for (std::size_t i = 0; i < a.nrows; i++) {
      std::ranges::fill(temp, -1);
      std::size_t ind1 = a.row_ptrs[i];
      std::size_t ind2 = a.row_ptrs[i + 1];
      for (std::size_t j = ind1; j < ind2; j++) {
        temp[a.cols[j]] = static_cast<int>(j);
      }
      for (std::size_t j = 0; j < bt.nrows; j++) {
        sum = 0;
        std::size_t ind3 = bt.row_ptrs[j];
        std::size_t ind4 = bt.row_ptrs[j + 1];
        for (std::size_t k = ind3; k < ind4; k++) {
          int aind = temp[bt.cols[k]];
          if (aind != -1) {
            sum += a.values[aind] * bt.values[k];
          }
        }

        if (std::fabs(sum) > kZERO) {
          result.cols.push_back(j);
          result.values.push_back(sum);
          nz++;
        }
      }
      result.row_ptrs.push_back(nz);
    }
    result.nz = nz;
    return result;
  }
};
}  // namespace guseva_crs
