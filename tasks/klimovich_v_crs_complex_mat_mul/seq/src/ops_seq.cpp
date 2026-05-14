#include "klimovich_v_crs_complex_mat_mul/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "klimovich_v_crs_complex_mat_mul/common/include/common.hpp"

namespace klimovich_v_crs_complex_mat_mul {
namespace {

void AccumulateRow(const CrsMatrix &lhs, const CrsMatrix &rhs, int row, std::vector<Cplx> &spa,
                   std::vector<int> &touched_by_row, std::vector<int> &touched_cols) {
  touched_cols.clear();
  for (int lp = lhs.row_offsets[row]; lp < lhs.row_offsets[row + 1]; ++lp) {
    const int k = lhs.col_indices[lp];
    const Cplx a_ik = lhs.data[lp];
    for (int rq = rhs.row_offsets[k]; rq < rhs.row_offsets[k + 1]; ++rq) {
      const int j = rhs.col_indices[rq];
      if (touched_by_row[j] != row) {
        touched_by_row[j] = row;
        touched_cols.push_back(j);
        spa[j] = a_ik * rhs.data[rq];
      } else {
        spa[j] += a_ik * rhs.data[rq];
      }
    }
  }
}

void FinalizeRow(CrsMatrix &result, std::vector<Cplx> &spa, std::vector<int> &touched_cols, int row) {
  std::ranges::sort(touched_cols);

  int kept = 0;
  for (const int j : touched_cols) {
    const Cplx v = spa[j];
    spa[j] = Cplx(0.0);
    if (std::abs(v.real()) > kZeroDropTol || std::abs(v.imag()) > kZeroDropTol) {
      result.col_indices.push_back(j);
      result.data.push_back(v);
      ++kept;
    }
  }

  result.row_offsets[row + 1] = result.row_offsets[row] + kept;
}

}  // namespace

KlimovichVCrsComplexMatMulSeq::KlimovichVCrsComplexMatMulSeq(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = CrsMatrix();
}

bool KlimovichVCrsComplexMatMulSeq::ValidationImpl() {
  const auto &lhs = std::get<0>(GetInput());
  const auto &rhs = std::get<1>(GetInput());
  return lhs.n_cols == rhs.n_rows;
}

bool KlimovichVCrsComplexMatMulSeq::PreProcessingImpl() {
  return true;
}

CrsMatrix KlimovichVCrsComplexMatMulSeq::MultiplyCrs(const CrsMatrix &lhs, const CrsMatrix &rhs) {
  CrsMatrix result(lhs.n_rows, rhs.n_cols);

  std::vector<Cplx> spa(static_cast<std::size_t>(rhs.n_cols), Cplx(0.0));
  std::vector<int> touched_by_row(static_cast<std::size_t>(rhs.n_cols), -1);
  std::vector<int> touched_cols;
  touched_cols.reserve(static_cast<std::size_t>(rhs.n_cols));

  for (int i = 0; i < lhs.n_rows; ++i) {
    AccumulateRow(lhs, rhs, i, spa, touched_by_row, touched_cols);
    FinalizeRow(result, spa, touched_cols, i);
  }

  return result;
}

bool KlimovichVCrsComplexMatMulSeq::RunImpl() {
  const auto &lhs = std::get<0>(GetInput());
  const auto &rhs = std::get<1>(GetInput());
  GetOutput() = MultiplyCrs(lhs, rhs);
  return true;
}

bool KlimovichVCrsComplexMatMulSeq::PostProcessingImpl() {
  return true;
}

}  // namespace klimovich_v_crs_complex_mat_mul
