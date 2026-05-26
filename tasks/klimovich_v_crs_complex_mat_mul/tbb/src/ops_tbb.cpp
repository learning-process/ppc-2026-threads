#include "klimovich_v_crs_complex_mat_mul/tbb/include/ops_tbb.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "klimovich_v_crs_complex_mat_mul/common/include/common.hpp"
#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/enumerable_thread_specific.h"
#include "oneapi/tbb/parallel_for.h"

namespace klimovich_v_crs_complex_mat_mul {
namespace {

struct RowStage {
  std::vector<int> cols;
  std::vector<Cplx> vals;
};

struct ThreadCtx {
  std::vector<Cplx> spa;
  std::vector<int> touched_by_row;
  std::vector<int> touched_cols;
};

void GustavsonRow(const CrsMatrix &lhs, const CrsMatrix &rhs, int row, ThreadCtx &ctx, RowStage &stage) {
  auto &spa = ctx.spa;
  auto &touched_by_row = ctx.touched_by_row;
  auto &touched_cols = ctx.touched_cols;
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

  std::ranges::sort(touched_cols);

  stage.cols.clear();
  stage.vals.clear();
  stage.cols.reserve(touched_cols.size());
  stage.vals.reserve(touched_cols.size());

  for (const int j : touched_cols) {
    const Cplx v = spa[j];
    spa[j] = Cplx(0.0, 0.0);
    if (std::abs(v.real()) > kZeroDropTol || std::abs(v.imag()) > kZeroDropTol) {
      stage.cols.push_back(j);
      stage.vals.push_back(v);
    }
  }
}

CrsMatrix Assemble(int rows, int cols, const std::vector<RowStage> &per_row) {
  CrsMatrix out(rows, cols);
  for (int i = 0; i < rows; ++i) {
    out.row_offsets[i + 1] = out.row_offsets[i] + static_cast<int>(per_row[i].cols.size());
  }
  out.col_indices.reserve(static_cast<std::size_t>(out.row_offsets[rows]));
  out.data.reserve(static_cast<std::size_t>(out.row_offsets[rows]));
  for (int i = 0; i < rows; ++i) {
    out.col_indices.insert(out.col_indices.end(), per_row[i].cols.begin(), per_row[i].cols.end());
    out.data.insert(out.data.end(), per_row[i].vals.begin(), per_row[i].vals.end());
  }
  return out;
}

}  // namespace

KlimovichVCrsComplexMatMulTbb::KlimovichVCrsComplexMatMulTbb(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = CrsMatrix();
}

bool KlimovichVCrsComplexMatMulTbb::ValidationImpl() {
  const auto &lhs = std::get<0>(GetInput());
  const auto &rhs = std::get<1>(GetInput());
  return lhs.n_cols == rhs.n_rows;
}

bool KlimovichVCrsComplexMatMulTbb::PreProcessingImpl() {
  return true;
}

CrsMatrix KlimovichVCrsComplexMatMulTbb::MultiplyCrs(const CrsMatrix &lhs, const CrsMatrix &rhs) {
  std::vector<RowStage> per_row(static_cast<std::size_t>(lhs.n_rows));

  oneapi::tbb::enumerable_thread_specific<ThreadCtx> tls([&rhs] {
    ThreadCtx c;
    c.spa.assign(static_cast<std::size_t>(rhs.n_cols), Cplx(0.0, 0.0));
    c.touched_by_row.assign(static_cast<std::size_t>(rhs.n_cols), -1);
    c.touched_cols.reserve(static_cast<std::size_t>(rhs.n_cols));
    return c;
  });

  oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<int>(0, lhs.n_rows), [&](const auto &range) {
    auto &ctx = tls.local();
    for (int i = range.begin(); i < range.end(); ++i) {
      GustavsonRow(lhs, rhs, i, ctx, per_row[i]);
    }
  });

  return Assemble(lhs.n_rows, rhs.n_cols, per_row);
}

bool KlimovichVCrsComplexMatMulTbb::RunImpl() {
  const auto &lhs = std::get<0>(GetInput());
  const auto &rhs = std::get<1>(GetInput());
  GetOutput() = MultiplyCrs(lhs, rhs);
  return true;
}

bool KlimovichVCrsComplexMatMulTbb::PostProcessingImpl() {
  return true;
}

}  // namespace klimovich_v_crs_complex_mat_mul
