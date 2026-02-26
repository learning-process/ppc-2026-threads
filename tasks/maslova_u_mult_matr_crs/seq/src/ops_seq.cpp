#include "maslova_u_mult_matr_crs/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

#include "maslova_u_mult_matr_crs/common/include/common.hpp"
#include "util/include/util.hpp"

namespace maslova_u_mult_matr_crs {

MaslovaUMultMatrSEQ::MaslovaUMultMatrSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool MaslovaUMultMatrSEQ::ValidationImpl() {
  const auto &A = std::get<0>(GetInput());
  const auto &B = std::get<1>(GetInput());
  if (A.cols != B.rows || A.rows <= 0 || B.cols <= 0) {
    return false;
  }
  if (A.row_ptr.size() != static_cast<size_t>(A.rows + 1)) {
    return false;
  }
  if (B.row_ptr.size() != static_cast<size_t>(B.rows + 1)) {
    return false;
  }
  return true;
}

bool MaslovaUMultMatrSEQ::PreProcessingImpl() {
  const auto &B = std::get<1>(GetInput());
  if (temp_row.size() != static_cast<size_t>(B.cols)) {
    temp_row.assign(B.cols, 0.0);
  }
  if (marker.size() != static_cast<size_t>(B.cols)) {
    marker.assign(B.cols, -1);
  }
  used_cols.reserve(B.cols);

  return true;
}

bool MaslovaUMultMatrSEQ::RunImpl() {
  const auto &A = std::get<0>(GetInput());
  const auto &B = std::get<1>(GetInput());
  auto &C = GetOutput();

  C.rows = A.rows;
  C.cols = B.cols;
  C.row_ptr.assign(A.rows + 1, 0);
  C.values.clear();
  C.col_ind.clear();

  C.values.reserve(A.values.size());
  C.col_ind.reserve(A.values.size());

  std::fill(marker.begin(), marker.end(), -1);
  used_cols.clear();

  for (int i = 0; i < A.rows; ++i) {
    for (int j = A.row_ptr[i]; j < A.row_ptr[i + 1]; ++j) {
      int col_A = A.col_ind[j];
      double val_A = A.values[j];

      for (int k = B.row_ptr[col_A]; k < B.row_ptr[col_A + 1]; ++k) {
        int col_B = B.col_ind[k];
        double val_B = B.values[k];

        if (marker[col_B] < i) {
          marker[col_B] = i;
          used_cols.push_back(col_B);
          temp_row[col_B] = val_A * val_B;
        } else {
          temp_row[col_B] += val_A * val_B;
        }
      }
    }

    if (!used_cols.empty()) {
      std::sort(used_cols.begin(), used_cols.end());

      for (int col_idx : used_cols) {
        double val = temp_row[col_idx];
        if (std::abs(val) > 1e-15) {
          C.values.push_back(val);
          C.col_ind.push_back(col_idx);
        }
        temp_row[col_idx] = 0.0;
      }
      used_cols.clear();
    }
    C.row_ptr[i + 1] = static_cast<int>(C.values.size());
  }

  return true;
}

bool MaslovaUMultMatrSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace maslova_u_mult_matr_crs
