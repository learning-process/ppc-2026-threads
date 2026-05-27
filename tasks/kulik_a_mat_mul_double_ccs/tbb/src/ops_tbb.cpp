#include "kulik_a_mat_mul_double_ccs/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <cstddef>
#include <limits>
#include <tuple>
#include <vector>

#include "kulik_a_mat_mul_double_ccs/common/include/common.hpp"

namespace kulik_a_mat_mul_double_ccs {

namespace {

inline void MatMultPhase1(size_t j, const CCS &a, const CCS &b, std::vector<size_t> &was,
                          std::vector<size_t> &col_nnz) {
  size_t count = 0;
  for (size_t k = b.col_ind[j]; k < b.col_ind[j + 1]; ++k) {
    const size_t b_row = b.row[k];
    for (size_t zc = a.col_ind[b_row]; zc < a.col_ind[b_row + 1]; ++zc) {
      const size_t a_row = a.row[zc];
      if (was[a_row] != j) {
        was[a_row] = j;
        ++count;
      }
    }
  }
  col_nnz[j] = count;
}

inline void MatMultPhase2(size_t j, const CCS &a, const CCS &b, CCS &c, size_t stamp, std::vector<size_t> &was,
                          std::vector<double> &accum, std::vector<size_t> &rows) {
  rows.clear();
  for (size_t k = b.col_ind[j]; k < b.col_ind[j + 1]; ++k) {
    const double b_val = b.value[k];
    const size_t b_row = b.row[k];
    for (size_t zc = a.col_ind[b_row]; zc < a.col_ind[b_row + 1]; ++zc) {
      const size_t a_row = a.row[zc];
      accum[a_row] += a.value[zc] * b_val;
      if (was[a_row] != stamp) {
        was[a_row] = stamp;
        rows.push_back(a_row);
      }
    }
  }
  std::ranges::sort(rows);
  size_t write = c.col_ind[j];
  for (const size_t i : rows) {
    c.row[write] = i;
    c.value[write] = accum[i];
    accum[i] = 0.0;
    ++write;
  }
}

}  // namespace

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

  std::vector<size_t> col_nnz(b.m, 0);

  tbb::enumerable_thread_specific<std::vector<size_t>> tls_was(
      [&]() { return std::vector<size_t>(a.n, std::numeric_limits<size_t>::max()); });

  tbb::parallel_for(tbb::blocked_range<size_t>(0, b.m), [&](const tbb::blocked_range<size_t> &r) {
    auto &was = tls_was.local();
    for (size_t j = r.begin(); j != r.end(); ++j) {
      MatMultPhase1(j, a, b, was, col_nnz);
    }
  });

  size_t total_nz = 0;
  for (size_t j = 0; j < b.m; ++j) {
    c.col_ind[j] = total_nz;
    total_nz += col_nnz[j];
  }
  c.col_ind[b.m] = total_nz;
  c.nz = total_nz;
  c.value.resize(total_nz);
  c.row.resize(total_nz);

  tbb::enumerable_thread_specific<std::vector<double>> tls_accum([&]() { return std::vector<double>(a.n, 0.0); });
  tbb::enumerable_thread_specific<std::vector<size_t>> tls_rows;

  tbb::parallel_for(tbb::blocked_range<size_t>(0, b.m), [&](const tbb::blocked_range<size_t> &r) {
    auto &was = tls_was.local();
    auto &accum = tls_accum.local();
    auto &rows = tls_rows.local();
    for (size_t j = r.begin(); j != r.end(); ++j) {
      MatMultPhase2(j, a, b, c, b.m + j, was, accum, rows);
    }
  });

  return true;
}

bool KulikAMatMulDoubleCcsTBB::PostProcessingImpl() {
  return true;
}

}  // namespace kulik_a_mat_mul_double_ccs
