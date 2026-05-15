#include "romanov_m_matrix_ccs/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "romanov_m_matrix_ccs/common/include/common.hpp"
#include "util/include/util.hpp"

namespace romanov_m_matrix_ccs {

RomanovMMatrixCCSOMP::RomanovMMatrixCCSOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool RomanovMMatrixCCSOMP::ValidationImpl() {
  return GetInput().first.cols_num == GetInput().second.rows_num;
}

bool RomanovMMatrixCCSOMP::PreProcessingImpl() {
  return true;
}

void RomanovMMatrixCCSOMP::MultiplyColumn(size_t col_index, const MatrixCCS &a, const MatrixCCS &b,
                                          std::vector<double> &temp_v, std::vector<size_t> &temp_r) {
  std::vector<double> accumulator(a.rows_num, 0.0);
  std::vector<bool> row_mask(a.rows_num, false);
  std::vector<size_t> active_rows;

  for (size_t kb = b.col_ptrs[col_index]; kb < b.col_ptrs[col_index + 1]; ++kb) {
    size_t k = b.row_inds[kb];
    double v_b = b.vals[kb];

    for (size_t ka = a.col_ptrs[k]; ka < a.col_ptrs[k + 1]; ++ka) {
      size_t i = a.row_inds[ka];
      if (!row_mask[i]) {
        row_mask[i] = true;
        active_rows.push_back(i);
      }
      accumulator[i] += a.vals[ka] * v_b;
    }
  }

  std::ranges::sort(active_rows);
  for (size_t row_idx : active_rows) {
    if (std::abs(accumulator[row_idx]) > 1e-12) {
      temp_v.push_back(accumulator[row_idx]);
      temp_r.push_back(row_idx);
    }
  }
}

bool RomanovMMatrixCCSOMP::RunImpl() {
  const auto &a = GetInput().first;
  const auto &b = GetInput().second;
  auto &c = GetOutput();

  c.rows_num = a.rows_num;
  c.cols_num = b.cols_num;
  c.col_ptrs.assign(c.cols_num + 1, 0);

  std::vector<std::vector<double>> temp_vals(c.cols_num);
  std::vector<std::vector<size_t>> temp_rows(c.cols_num);

#pragma omp parallel num_threads(ppc::util::GetNumThreads()) default(none) shared(a, b, temp_vals, temp_rows)
  {
#pragma omp for schedule(dynamic)
    for (size_t j = 0; j < b.cols_num; ++j) {
      MultiplyColumn(j, a, b, temp_vals[j], temp_rows[j]);
    }
  }

  size_t total_nnz = 0;
  for (size_t j = 0; j < b.cols_num; ++j) {
    c.col_ptrs[j] = total_nnz;
    total_nnz += temp_vals[j].size();
  }
  c.col_ptrs[b.cols_num] = total_nnz;
  c.nnz = total_nnz;

  c.vals.reserve(total_nnz);
  c.row_inds.reserve(total_nnz);
  for (size_t j = 0; j < b.cols_num; ++j) {
    c.vals.insert(c.vals.end(), temp_vals[j].begin(), temp_vals[j].end());
    c.row_inds.insert(c.row_inds.end(), temp_rows[j].begin(), temp_rows[j].end());
  }

  return true;
}

bool RomanovMMatrixCCSOMP::PostProcessingImpl() {
  return true;
}

}  // namespace romanov_m_matrix_ccs
