#include "zagryadskov_m_complex_spmm_ccs/omp/include/ops_omp.hpp"

#include <omp.h>

#include <cmath>
#include <complex>
#include <vector>

#include "util/include/util.hpp"
#include "zagryadskov_m_complex_spmm_ccs/common/include/common.hpp"

namespace zagryadskov_m_complex_spmm_ccs {

ZagryadskovMComplexSpMMCCSOMP::ZagryadskovMComplexSpMMCCSOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = CCS();
}

/*
void ZagryadskovMComplexSpMMCCSOMP::SpMMkernel(const CCS &a, const CCS &b, const std::complex<double> &zero, double eps,
                                               int num_threads, std::vector<std::vector<int>> &t_row_ind,
                                               std::vector<std::vector<std::complex<double>>> &t_values,
                                               std::vector<std::vector<int>> &t_col_ptr) {
  int tid = omp_get_thread_num();
  int jstart = (tid * b.n) / num_threads;
  int jend = ((tid + 1) * b.n) / num_threads;

  t_col_ptr[tid].assign(jend - jstart + 1, 0);

  std::vector<int> rows;
  std::vector<int> marker(a.m, -1);
  std::vector<std::complex<double>> acc(a.m, zero);

  for (int j = jstart; j < jend; ++j) {
    rows.clear();

    for (int k = b.col_ptr[j]; k < b.col_ptr[j + 1]; ++k) {
      std::complex<double> tmpval = b.values[k];
      int btmpind = b.row_ind[k];

      for (int zp = a.col_ptr[btmpind]; zp < a.col_ptr[btmpind + 1]; ++zp) {
        int atmpind = a.row_ind[zp];
        acc[atmpind] += tmpval * a.values[zp];
        if (marker[atmpind] != j) {
          rows.push_back(atmpind);
          marker[atmpind] = j;
        }
      }
    }

    for (int tmpind : rows) {
      if (std::abs(acc[tmpind]) > eps) {
        t_values[tid].push_back(acc[tmpind]);
        t_row_ind[tid].push_back(tmpind);
        ++t_col_ptr[tid][j - jstart + 1];
      }
      acc[tmpind] = zero;
    }
  }
}

void ZagryadskovMComplexSpMMCCSOMP::SpMM(const CCS &a, const CCS &b, CCS &c) {
  c.m = a.m;
  c.n = b.n;
  c.col_ptr.assign(b.n + 1, 0);
  c.row_ind.clear();
  c.values.clear();
  std::complex<double> zero(0.0, 0.0);
  const double eps = 1e-14;
  const int num_threads = ppc::util::GetNumThreads();

  std::vector<std::vector<int>> t_row_ind(num_threads);
  std::vector<std::vector<std::complex<double>>> t_values(num_threads);
  std::vector<std::vector<int>> t_col_ptr(num_threads);

#pragma omp parallel default(none) shared(zero, eps, num_threads, t_row_ind, t_values, t_col_ptr, a, b) \
    num_threads(ppc::util::GetNumThreads())
  {
    SpMMkernel(a, b, zero, eps, num_threads, t_row_ind, t_values, t_col_ptr);
  }
  for (int tid = 0; tid < num_threads; ++tid) {
    int jstart = (tid * b.n) / num_threads;
    int jend = ((tid + 1) * b.n) / num_threads;

    for (int j = jstart; j < jend; ++j) {
      c.col_ptr[j + 1] = c.col_ptr[j] + t_col_ptr[tid][j - jstart + 1];
    }

    c.row_ind.insert(c.row_ind.end(), t_row_ind[tid].begin(), t_row_ind[tid].end());
    c.values.insert(c.values.end(), t_values[tid].begin(), t_values[tid].end());
  }
}
  */

void ZagryadskovMComplexSpMMCCSOMP::SpMMkernel(const CCS &a, const CCS &b, const std::complex<double> &zero, double eps,
                                               std::vector<std::vector<int>> &col_rows,
                                               std::vector<std::vector<std::complex<double>>> &col_vals) {
#pragma omp parallel num_threads(ppc::util::GetNumThreads())
  {
    std::vector<std::complex<double>> acc(a.m, zero);
    std::vector<int> marker(a.m, -1);
    std::vector<int> rows;
    rows.reserve(a.m);

#pragma omp for schedule(static)
    for (int j = 0; j < b.n; ++j) {
      rows.clear();

      for (int k = b.col_ptr[j]; k < b.col_ptr[j + 1]; ++k) {
        std::complex<double> bval = b.values[k];
        int bcol = b.row_ind[k];

        for (int zp = a.col_ptr[bcol]; zp < a.col_ptr[bcol + 1]; ++zp) {
          int arow = a.row_ind[zp];
          acc[arow] += bval * a.values[zp];
          if (marker[arow] != j) {
            rows.push_back(arow);
            marker[arow] = j;
          }
        }
      }

      for (int r : rows) {
        if (std::abs(acc[r]) > eps) {
          col_rows[j].push_back(r);
          col_vals[j].push_back(acc[r]);
        }
        acc[r] = zero;
      }
    }
  }
}

void ZagryadskovMComplexSpMMCCSOMP::SpMM(const CCS &a, const CCS &b, CCS &c) {
  c.m = a.m;
  c.n = b.n;
  const double eps = 1e-14;
  const std::complex<double> zero(0.0, 0.0);
  std::vector<std::vector<int>> col_rows(b.n);
  std::vector<std::vector<std::complex<double>>> col_vals(b.n);

  SpMMkernel(a, b, zero, eps, col_rows, col_vals);

  c.col_ptr.resize(b.n + 1);
  c.col_ptr[0] = 0;
  for (int j = 0; j < b.n; ++j) {
    c.col_ptr[j + 1] = c.col_ptr[j] + static_cast<int>(col_rows[j].size());
  }

  int nnz = c.col_ptr[b.n];
  c.row_ind.resize(nnz);
  c.values.resize(nnz);

  for (int j = 0; j < b.n; ++j) {
    int offset = c.col_ptr[j];
    std::copy(col_rows[j].begin(), col_rows[j].end(), c.row_ind.begin() + offset);
    std::copy(col_vals[j].begin(), col_vals[j].end(), c.values.begin() + offset);
  }
}

bool ZagryadskovMComplexSpMMCCSOMP::ValidationImpl() {
  const CCS &a = std::get<0>(GetInput());
  const CCS &b = std::get<1>(GetInput());
  return a.n == b.m;
}

bool ZagryadskovMComplexSpMMCCSOMP::PreProcessingImpl() {
  return true;
}

bool ZagryadskovMComplexSpMMCCSOMP::RunImpl() {
  const CCS &a = std::get<0>(GetInput());
  const CCS &b = std::get<1>(GetInput());
  CCS &c = GetOutput();

  ZagryadskovMComplexSpMMCCSOMP::SpMM(a, b, c);

  return true;
}

bool ZagryadskovMComplexSpMMCCSOMP::PostProcessingImpl() {
  return true;
}

}  // namespace zagryadskov_m_complex_spmm_ccs
