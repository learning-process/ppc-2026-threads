#include "../include/ops_tbb.hpp"

#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>

#include "../../common/include/common.hpp"
#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/parallel_for.h"

namespace sabutay_sparse_complex_ccs_mult_tbb {

namespace {

bool IsValidCCS(const CCS &matrix) {
  if (matrix.m < 0 || matrix.n < 0) {
    return false;
  }
  if (matrix.col_ptr.size() != static_cast<std::size_t>(matrix.n + 1)) {
    return false;
  }
  if (matrix.row_ind.size() != matrix.values.size()) {
    return false;
  }
  if (matrix.col_ptr.empty() || matrix.col_ptr[0] != 0) {
    return false;
  }
  if (matrix.col_ptr.back() != static_cast<int>(matrix.row_ind.size())) {
    return false;
  }
  for (int j = 0; j < matrix.n; ++j) {
    if (matrix.col_ptr[j] > matrix.col_ptr[j + 1]) {
      return false;
    }
  }
  for (int row : matrix.row_ind) {
    if (row < 0 || row >= matrix.m) {
      return false;
    }
  }
  return true;
}

}  // namespace

SabutayASparseComplexCcsMultTBB::SabutayASparseComplexCcsMultTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = CCS();
}

void SabutayASparseComplexCcsMultTBB::SpMM(const CCS &a, const CCS &b, CCS &c) {
  c.m = a.m;
  c.n = b.n;
  c.col_ptr.assign(b.n + 1, 0);
  c.row_ind.clear();
  c.values.clear();

  std::vector<std::vector<int>> local_row_ind(b.n);
  std::vector<std::vector<std::complex<double>>> local_values(b.n);
  std::vector<int> local_sizes(b.n, 0);
  const double eps = 1e-14;

  oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<int>(0, b.n), [&](const oneapi::tbb::blocked_range<int> &range) {
    std::complex<double> zero(0.0, 0.0);
    std::vector<int> rows;
    std::vector<int> marker(a.m, -1);
    std::vector<std::complex<double>> acc(a.m);

    for (int j = range.begin(); j < range.end(); ++j) {
      rows.clear();
      local_row_ind[j].clear();
      local_values[j].clear();

      for (int k = b.col_ptr[j]; k < b.col_ptr[j + 1]; ++k) {
        std::complex<double> tmpval = b.values[k];
        int btmpind = b.row_ind[k];

        for (int zp = a.col_ptr[btmpind]; zp < a.col_ptr[btmpind + 1]; ++zp) {
          int atmpind = a.row_ind[zp];
          acc[atmpind] += tmpval * a.values[zp];
          if (marker[atmpind] == -1) {
            rows.push_back(atmpind);
            marker[atmpind] = 1;
          }
        }
      }

      local_row_ind[j].reserve(rows.size());
      local_values[j].reserve(rows.size());

      for (int tmpind : rows) {
        if (std::abs(acc[tmpind]) > eps) {
          local_values[j].push_back(acc[tmpind]);
          local_row_ind[j].push_back(tmpind);
        }
        acc[tmpind] = zero;
        marker[tmpind] = -1;
      }

      local_sizes[j] = static_cast<int>(local_values[j].size());
    }
  });

  for (int j = 0; j < b.n; ++j) {
    c.col_ptr[j + 1] = c.col_ptr[j] + local_sizes[j];
  }

  c.row_ind.reserve(c.col_ptr.back());
  c.values.reserve(c.col_ptr.back());

  for (int j = 0; j < b.n; ++j) {
    c.row_ind.insert(c.row_ind.end(), local_row_ind[j].begin(), local_row_ind[j].end());
    c.values.insert(c.values.end(), local_values[j].begin(), local_values[j].end());
  }
}

bool SabutayASparseComplexCcsMultTBB::ValidationImpl() {
  const CCS &a = std::get<0>(GetInput());
  const CCS &b = std::get<1>(GetInput());
  return IsValidCCS(a) && IsValidCCS(b) && a.n == b.m;
}

bool SabutayASparseComplexCcsMultTBB::PreProcessingImpl() {
  return true;
}

bool SabutayASparseComplexCcsMultTBB::RunImpl() {
  const CCS &a = std::get<0>(GetInput());
  const CCS &b = std::get<1>(GetInput());
  CCS &c = GetOutput();

  SpMM(a, b, c);

  return true;
}

bool SabutayASparseComplexCcsMultTBB::PostProcessingImpl() {
  return true;
}

}  // namespace sabutay_sparse_complex_ccs_mult_tbb
