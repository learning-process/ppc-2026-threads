#pragma once
#include <cmath>
#include <cstddef>
#include <functional>
#include <numeric>
#include <thread>
#include <vector>

#include "guseva_crs/common/include/common.hpp"
#include "guseva_crs/common/include/multiplier.hpp"
#include "util/include/util.hpp"

namespace guseva_crs {

class MultiplierStl : public Multiplier {
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

  static void ProcessRows(std::size_t i, const CRS &a, const CRS &bt, std::vector<std::vector<std::size_t>> &columns,
                          std::vector<std::vector<double>> &values, std::vector<std::size_t> &row_index) {
    std::size_t n = a.nrows;
    std::vector<int> temp(n);

    for (int &l : temp) {
      l = -1;
    }

    std::size_t ind1 = a.row_ptrs[i];
    std::size_t ind2 = a.row_ptrs[i + 1];
    for (std::size_t j = ind1; j < ind2; j++) {
      std::size_t col = a.cols[j];
      temp[col] = static_cast<int>(j);
    }

    for (std::size_t j = 0; j < n; j++) {
      double sum = 0;
      std::size_t ind3 = bt.row_ptrs[j];
      std::size_t ind4 = bt.row_ptrs[j + 1];

      PerformCalculation(0, ind3, ind4, a, bt, sum, temp);

      if (std::fabs(sum) > kZERO) {
        columns[i].push_back(j);
        values[i].push_back(sum);
        row_index[i]++;
      }
    }
  }

  static void ProcessRowsRange(const std::vector<std::size_t> &indices, std::size_t start, std::size_t end,
                               const CRS &a, const CRS &bt, std::vector<std::vector<std::size_t>> &columns,
                               std::vector<std::vector<double>> &values, std::vector<std::size_t> &row_index) {
    for (std::size_t idx = start; idx < end; ++idx) {
      std::size_t i = indices[idx];
      ProcessRows(i, a, bt, columns, values, row_index);
    }
  }

 public:
  [[nodiscard]] CRS Multiply(const CRS &a, const CRS &b) const override {
    std::size_t n = a.nrows;

    auto bt = this->Transpose(b);

    std::vector<std::vector<std::size_t>> columns(n);
    std::vector<std::vector<double>> values(n);
    std::vector<std::size_t> row_index(n + 1, 0);

    std::vector<std::size_t> indices(n);
#ifdef __APPLE__
    std::iota(indices.begin(), indices.end(), 0);
#else
    std::ranges::iota(indices, 0);
#endif

    std::size_t num_threads = ppc::util::GetNumThreads();
    if (num_threads == 0) {
      num_threads = 2;
    }

    std::vector<std::thread> threads;
    std::size_t chunk_size = n / num_threads;
    std::size_t remainder = n % num_threads;

    std::size_t start = 0;
    for (std::size_t thread = 0; thread < num_threads; ++thread) {
      std::size_t end = start + chunk_size + (thread < remainder ? 1 : 0);

      threads.emplace_back(ProcessRowsRange, std::ref(indices), start, end, std::cref(a), std::cref(bt),
                           std::ref(columns), std::ref(values), std::ref(row_index));
      start = end;
    }

    for (auto &thread : threads) {
      thread.join();
    }

    std::size_t nz = 0;
    for (std::size_t i = 0; i < n; i++) {
      std::size_t tmp = row_index[i];
      row_index[i] = nz;
      nz += tmp;
    }
    row_index[n] = nz;

    CRS result;
    result.row_ptrs = row_index;
    result.nrows = n;
    result.ncols = n;

    for (std::size_t i = 0; i < n; i++) {
      result.cols.insert(result.cols.end(), columns[i].begin(), columns[i].end());
      result.values.insert(result.values.end(), values[i].begin(), values[i].end());
    }

    result.nz = result.values.size();
    return result;
  }
};

}  // namespace guseva_crs
