#include "zagryadskov_m_complex_spmm_ccs/stl/include/ops_stl.hpp"

#include <complex>
#include <functional>
#include <thread>
#include <vector>

#include "util/include/util.hpp"
#include "zagryadskov_m_complex_spmm_ccs/common/include/common.hpp"

namespace zagryadskov_m_complex_spmm_ccs {

ZagryadskovMComplexSpMMCCSSTL::ZagryadskovMComplexSpMMCCSSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = CCS();
}

void ZagryadskovMComplexSpMMCCSSTL::SpMMSymbolic(const CCS &a, const CCS &b, std::vector<int> &col_ptr, int jstart,
                                                 int jend) {
  std::vector<int> marker(a.m, -1);

  for (int j = jstart; j < jend; ++j) {
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
    col_ptr[j + 1] += count;
  }
}

void ZagryadskovMComplexSpMMCCSSTL::SpMMKernel(const CCS &a, const CCS &b, CCS &c, const std::complex<double> &zero,
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

void ZagryadskovMComplexSpMMCCSSTL::SpMMNumeric(const CCS &a, const CCS &b, CCS &c, const std::complex<double> &zero,
                                                double eps, int jstart, int jend) {
  std::vector<int> marker(a.m, -1);
  std::vector<std::complex<double>> acc(a.m, zero);
  std::vector<int> rows;

  for (int j = jstart; j < jend; ++j) {
    SpMMKernel(a, b, c, zero, eps, rows, acc, marker, j);
  }
}

void ZagryadskovMComplexSpMMCCSSTL::SpMM(const CCS &a, const CCS &b, CCS &c) {
  c.m = a.m;
  c.n = b.n;
  const int num_threads = ppc::util::GetNumThreads();
  std::vector<std::thread> threads(num_threads);

  std::complex<double> zero(0.0, 0.0);
  const double eps = 1e-14;
  c.col_ptr.assign(c.n + 1, 0);

  for (int tid = 0; tid < num_threads; ++tid) {
    int jstart = (tid * b.n) / num_threads;
    int jend = ((tid + 1) * b.n) / num_threads;
    threads[tid] = std::thread(SpMMSymbolic, std::cref(a), std::cref(b), std::ref(c.col_ptr), jstart, jend);
  }
  for (auto &th : threads) {
    th.join();
  }

  for (int j = 0; j < c.n; ++j) {
    c.col_ptr[j + 1] += c.col_ptr[j];
  }
  int nnz = c.col_ptr[b.n];
  c.row_ind.resize(nnz);
  c.values.resize(nnz);

  for (int tid = 0; tid < num_threads; ++tid) {
    int jstart = (tid * b.n) / num_threads;
    int jend = ((tid + 1) * b.n) / num_threads;
    threads[tid] =
        std::thread(SpMMNumeric, std::cref(a), std::cref(b), std::ref(c), std::cref(zero), eps, jstart, jend);
  }
  for (auto &th : threads) {
    th.join();
  }
}

bool ZagryadskovMComplexSpMMCCSSTL::ValidationImpl() {
  const CCS &a = std::get<0>(GetInput());
  const CCS &b = std::get<1>(GetInput());
  return a.n == b.m;
}

bool ZagryadskovMComplexSpMMCCSSTL::PreProcessingImpl() {
  return true;
}

bool ZagryadskovMComplexSpMMCCSSTL::RunImpl() {
  const CCS &a = std::get<0>(GetInput());
  const CCS &b = std::get<1>(GetInput());
  CCS &c = GetOutput();

  SpMM(a, b, c);

  return true;
}

bool ZagryadskovMComplexSpMMCCSSTL::PostProcessingImpl() {
  return true;
}

}  // namespace zagryadskov_m_complex_spmm_ccs
