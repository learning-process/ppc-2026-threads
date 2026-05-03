#include "kulik_a_mat_mul_double_ccs/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/tbb.h>

#include <algorithm>
#include <cstddef>
#include <tuple>
#include <vector>

#include "kulik_a_mat_mul_double_ccs/common/include/common.hpp"
#include "tbb/parallel_for.h"

namespace kulik_a_mat_mul_double_ccs {

namespace {

inline void ProcessColumn(size_t j, const CCS &a, const CCS &b, std::vector<double> &accum,
                          std::vector<bool> &nz_elem_rows, std::vector<size_t> &nnz_rows,
                          std::vector<std::vector<double>> &local_values,
                          std::vector<std::vector<size_t>> &local_rows) {
  for (size_t k = b.col_ind[j]; k < b.col_ind[j + 1]; ++k) {
    size_t ind = b.row[k];
    double b_val = b.value[k];
    for (size_t zc = a.col_ind[ind]; zc < a.col_ind[ind + 1]; ++zc) {
      size_t i = a.row[zc];
      double a_val = a.value[zc];

      accum[i] += a_val * b_val;
      if (!nz_elem_rows[i]) {
        nz_elem_rows[i] = true;
        nnz_rows.push_back(i);
      }
    }
  }

  std::ranges::sort(nnz_rows);

  for (size_t i : nnz_rows) {
    if (accum[i] != 0.0) {
      local_rows[j].push_back(i);
      local_values[j].push_back(accum[i]);
    }
    accum[i] = 0.0;
    nz_elem_rows[i] = false;
  }
  nnz_rows.clear();
}

inline void CopyColumn(size_t j, CCS &c, const std::vector<std::vector<double>> &local_values,
                       const std::vector<std::vector<size_t>> &local_rows) {
  size_t offset = c.col_ind[j];
  size_t col_nz = local_values[j].size();
  for (size_t k = 0; k < col_nz; ++k) {
    c.value[offset + k] = local_values[j][k];
    c.row[offset + k] = local_rows[j][k];
  }
}

}  // namespace

struct ThreadLocalData {
  std::vector<double> accum;
  std::vector<bool> nz_elem_rows;
  std::vector<size_t> nnz_rows;
};

KulikAMatMulDoubleCcsTBB::KulikAMatMulDoubleCcsTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool KulikAMatMulDoubleCcsTBB::ValidationImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());
  return (a.m == b.n);
}

bool KulikAMatMulDoubleCcsTBB::PreProcessingImpl() {
  return true;
}

bool KulikAMatMulDoubleCcsTBB::RunImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());
  OutType &c = GetOutput();

  c.n = a.n;
  c.m = b.m;
  c.col_ind.assign(c.m + 1, 0);

  std::vector<std::vector<double>> local_values(b.m);
  std::vector<std::vector<size_t>> local_rows(b.m);

  ThreadLocalData exemplar;
  exemplar.accum.assign(a.n, 0.0);
  exemplar.nz_elem_rows.assign(a.n, false);

  tbb::enumerable_thread_specific<ThreadLocalData> tls(exemplar);

  tbb::parallel_for(tbb::blocked_range<size_t>(0, b.m), [&](const tbb::blocked_range<size_t> &r) {
    auto &t_data = tls.local();

    for (size_t j = r.begin(); j != r.end(); ++j) {
      ProcessColumn(j, a, b, t_data.accum, t_data.nz_elem_rows, t_data.nnz_rows, local_values, local_rows);
    }
  });

  size_t total_nz = 0;
  for (size_t j = 0; j < b.m; ++j) {
    c.col_ind[j] = total_nz;
    total_nz += local_values[j].size();
  }
  c.col_ind[b.m] = total_nz;
  c.nz = total_nz;

  c.value.resize(total_nz);
  c.row.resize(total_nz);

  tbb::parallel_for(tbb::blocked_range<size_t>(0, b.m), [&](const tbb::blocked_range<size_t> &r) {
    for (size_t j = r.begin(); j != r.end(); ++j) {
      CopyColumn(j, c, local_values, local_rows);
    }
  });

  return true;
}

bool KulikAMatMulDoubleCcsTBB::PostProcessingImpl() {
  return true;
}

}  // namespace kulik_a_mat_mul_double_ccs
