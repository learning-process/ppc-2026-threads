#include "zagryadskov_m_complex_spmm_ccs/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/tbb.h>

#include <complex>
#include <vector>

#include "tbb/parallel_for.h"
#include "zagryadskov_m_complex_spmm_ccs/common/include/common.hpp"

namespace zagryadskov_m_complex_spmm_ccs {

ZagryadskovMComplexSpMMCCSTBB::ZagryadskovMComplexSpMMCCSTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = CCS();
}

void ZagryadskovMComplexSpMMCCSTBB::SpMMSymbolic(const CCS &a, const CCS &b, std::vector<int> &col_ptr) {
  const int m = a.m;
  const int n = b.n;

  std::vector<int> col_counts(n, 0);

  tbb::enumerable_thread_specific<std::vector<int>> tls_marker([&]() { return std::vector<int>(m, -1); });

  tbb::parallel_for(tbb::blocked_range<int>(0, n, 64), [&](const tbb::blocked_range<int> &r) {
    auto &marker = tls_marker.local();

    for (int j = r.begin(); j < r.end(); ++j) {
      int count = 0;

      for (int k = b.col_ptr[j]; k < b.col_ptr[j + 1]; ++k) {
        int b_row = b.row_ind[k];

        for (int zp = a.col_ptr[b_row]; zp < a.col_ptr[b_row + 1]; ++zp) {
          int a_row = a.row_ind[zp];

          if (marker[a_row] != j) {
            marker[a_row] = j;
            ++count;
          }
        }
      }

      col_counts[j] = count;
    }
  });
  col_ptr.resize(n + 1);
  col_ptr[0] = 0;
  for (int j = 0; j < n; ++j) {
    col_ptr[j + 1] = col_ptr[j] + col_counts[j];
  }
}

void ZagryadskovMComplexSpMMCCSTBB::SpMMKernel(const CCS &a, const CCS &b, CCS &c, const std::complex<double> &zero,
                                               double eps, std::vector<int> &rows,
                                               std::vector<std::complex<double>> &acc, std::vector<int> &marker,
                                               int j) {
  rows.clear();

  int write_ptr = c.col_ptr[j];

  for (int k = b.col_ptr[j]; k < b.col_ptr[j + 1]; ++k) {
    std::complex<double> tmpval = b.values[k];
    int b_row = b.row_ind[k];

    for (int zp = a.col_ptr[b_row]; zp < a.col_ptr[b_row + 1]; ++zp) {
      int a_row = a.row_ind[zp];

      acc[a_row] += tmpval * a.values[zp];

      if (marker[a_row] != j) {
        marker[a_row] = j;
        rows.push_back(a_row);
      }
    }
  }

  for (int r_idx : rows) {
    if (std::norm(acc[r_idx]) > eps * eps) {
      c.row_ind[write_ptr] = r_idx;
      c.values[write_ptr] = acc[r_idx];
      ++write_ptr;
    }
    acc[r_idx] = zero;
  }
}

void ZagryadskovMComplexSpMMCCSTBB::SpMMNumeric(const CCS &a, const CCS &b, CCS &c, const std::complex<double> &zero,
                                                double eps) {
  const int m = a.m;
  const int n = b.n;

  tbb::enumerable_thread_specific<std::vector<int>> tls_marker([&]() { return std::vector<int>(m, -1); });

  tbb::enumerable_thread_specific<std::vector<std::complex<double>>> tls_acc(
      [&]() { return std::vector<std::complex<double>>(m, zero); });

  tbb::parallel_for(tbb::blocked_range<int>(0, n, 64), [&](const tbb::blocked_range<int> &r) {
    auto &marker = tls_marker.local();
    auto &acc = tls_acc.local();

    std::vector<int> rows;

    for (int j = r.begin(); j < r.end(); ++j) {
      SpMMKernel(a, b, c, zero, eps, rows, acc, marker, j);
    }
  });
}

void ZagryadskovMComplexSpMMCCSTBB::SpMM(const CCS &a, const CCS &b, CCS &c) {
  c.m = a.m;
  c.n = b.n;

  std::complex<double> zero(0.0, 0.0);
  const double eps = 1e-14;

  SpMMSymbolic(a, b, c.col_ptr);

  int nnz = c.col_ptr[b.n];
  c.row_ind.resize(nnz);
  c.values.resize(nnz);

  SpMMNumeric(a, b, c, zero, eps);
}

bool ZagryadskovMComplexSpMMCCSTBB::ValidationImpl() {
  const CCS &a = std::get<0>(GetInput());
  const CCS &b = std::get<1>(GetInput());
  return a.n == b.m;
}

bool ZagryadskovMComplexSpMMCCSTBB::PreProcessingImpl() {
  return true;
}

bool ZagryadskovMComplexSpMMCCSTBB::RunImpl() {
  const CCS &a = std::get<0>(GetInput());
  const CCS &b = std::get<1>(GetInput());
  CCS &c = GetOutput();

  ZagryadskovMComplexSpMMCCSTBB::SpMM(a, b, c);

  return true;
}

bool ZagryadskovMComplexSpMMCCSTBB::PostProcessingImpl() {
  return true;
}

}  // namespace zagryadskov_m_complex_spmm_ccs
