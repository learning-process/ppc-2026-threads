#pragma once
#include <omp.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include "guseva_crs/common/include/common.hpp"
#include "guseva_crs/common/include/multiplier.hpp"

namespace guseva_crs {

class MultiplierOmp : public Multiplier {
  static void PerformCalculation(std::size_t k, std::size_t ind3, std::size_t ind4, const CRS &a, const CRS &bt,
                                 double &sum, std::vector<int> &temp) {
    for (k = ind3; k < ind4; k++) {
      std::size_t bcol = bt.cols[k];
      int aind = temp[bcol];
      if (aind != -1) {
        sum += a.values[aind] * bt.values[k];
      }
    }
  }

 public:
  [[nodiscard]] CRS Multiply(const CRS &a, const CRS &b) const override {
    std::size_t n = a.nrows;
    std::size_t i = 0;

    auto bt = this->Transpose(b);

    std::vector<std::vector<std::size_t>> columns(n);
    std::vector<std::vector<double>> values(n);
    std::vector<std::size_t> row_index(n + 1, 0);

#pragma omp parallel shared(n, a, bt, columns, values, row_index) private(i) default(none)
    {
      std::size_t j = 0;
      std::size_t k = 0;
      std::vector<int> temp(n);
#pragma omp for private(j, k)
      for (i = 0; i < n; i++) {
        for (int &l : temp) {
          l = -1;
        }
        std::size_t ind1 = a.row_ptrs[i];
        std::size_t ind2 = a.row_ptrs[i + 1];
        for (j = ind1; j < ind2; j++) {
          std::size_t col = a.cols[j];
          temp[col] = static_cast<int>(j);
        }
        for (j = 0; j < n; j++) {
          double sum = 0;
          std::size_t ind3 = bt.row_ptrs[j];
          std::size_t ind4 = bt.row_ptrs[j + 1];
          k = 0;
          PerformCalculation(k, ind3, ind4, a, bt, sum, temp);
          if (std::fabs(sum) > kZERO) {
            columns[i].push_back(j);
            values[i].push_back(sum);
            row_index[i]++;
          }
        }
      }
      temp.clear();
    }

    std::size_t nz = 0;
    for (i = 0; i < n; i++) {
      std::size_t tmp = row_index[i];
      row_index[i] = nz;
      nz += tmp;
    }

    CRS result;
    for (i = 0; i < n; i++) {
      result.cols.insert(result.cols.end(), columns[i].begin(), columns[i].end());
      result.values.insert(result.values.end(), values[i].begin(), values[i].end());
    }
    result.row_ptrs = row_index;
    result.nz = result.values.size();
    result.ncols = n;
    result.nrows = n;

    return result;
  }
};
}  // namespace guseva_crs
