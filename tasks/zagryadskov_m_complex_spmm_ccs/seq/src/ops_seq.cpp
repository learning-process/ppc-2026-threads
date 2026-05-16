#include "zagryadskov_m_complex_spmm_ccs/seq/include/ops_seq.hpp"

#include <cmath>
#include <complex>
#include <tuple>
#include <vector>

#include "zagryadskov_m_complex_spmm_ccs/common/include/common.hpp"

namespace zagryadskov_m_complex_spmm_ccs {

ZagryadskovMComplexSpMMCCSSEQ::ZagryadskovMComplexSpMMCCSSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = CCS();
}

void ZagryadskovMComplexSpMMCCSSEQ::SpMM(const CCS &a, const CCS &b, CCS &c) {
  c.m = a.m;
  c.n = b.n;
  c.col_ptr.assign(b.n + 1, 0);
  c.row_ind.clear();
  c.values.clear();
  std::complex<double> zero(0.0, 0.0);
  std::vector<int> rows;
  std::vector<int> marker(a.m, -1);
  std::vector<std::complex<double>> acc(a.m);
  const double eps = 1e-14;

  for (int j = 0; j < b.n; ++j) {
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
        c.values.push_back(acc[tmpind]);
        c.row_ind.push_back(tmpind);
        ++c.col_ptr[j + 1];
      }
      acc[tmpind] = zero;
    }

    c.col_ptr[j + 1] += c.col_ptr[j];
  }
}

bool ZagryadskovMComplexSpMMCCSSEQ::ValidationImpl() {
  const CCS &a = std::get<0>(GetInput());
  const CCS &b = std::get<1>(GetInput());
  return a.n == b.m;
}

bool ZagryadskovMComplexSpMMCCSSEQ::PreProcessingImpl() {
  return true;
}

bool ZagryadskovMComplexSpMMCCSSEQ::RunImpl() {
  const CCS &a = std::get<0>(GetInput());
  const CCS &b = std::get<1>(GetInput());
  CCS &c = GetOutput();

  SpMM(a, b, c);

  return true;
}

bool ZagryadskovMComplexSpMMCCSSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace zagryadskov_m_complex_spmm_ccs
