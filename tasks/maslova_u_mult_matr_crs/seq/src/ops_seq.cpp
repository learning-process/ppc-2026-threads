#include "maslova_u_mult_matr_crs/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "maslova_u_mult_matr_crs/common/include/common.hpp"

namespace maslova_u_mult_matr_crs {

MaslovaUMultMatrSEQ::MaslovaUMultMatrSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool MaslovaUMultMatrSEQ::ValidationImpl() {
  const auto &input = GetInput();
  const auto &a = std::get<0>(input);
  const auto &b = std::get<1>(input);

  if (a.cols != b.rows || a.rows <= 0 || b.cols <= 0) {
    return false;
  }
  if (a.row_ptr.size() != static_cast<size_t>(a.rows) + 1) {
    return false;
  }
  if (b.row_ptr.size() != static_cast<size_t>(b.rows) + 1) {
    return false;
  }
  return true;
}

bool MaslovaUMultMatrSEQ::PreProcessingImpl() {
  const auto &b = std::get<1>(GetInput());
  if (temp_row_.size() != static_cast<size_t>(b.cols)) {
    temp_row_.assign(static_cast<size_t>(b.cols), 0.0);
  }
  if (marker_.size() != static_cast<size_t>(b.cols)) {
    marker_.assign(static_cast<size_t>(b.cols), -1);
  }
  used_cols_.reserve(static_cast<size_t>(b.cols));

  return true;
}

bool MaslovaUMultMatrSEQ::RunImpl() {
  const auto &input = GetInput();
  const auto &a = std::get<0>(input);
  const auto &b = std::get<1>(input);
  auto &c = GetOutput();

  c.rows = a.rows;
  c.cols = b.cols;
  c.row_ptr.assign(static_cast<size_t>(a.rows) + 1, 0);
  c.values.clear();
  c.col_ind.clear();

  c.values.reserve(a.values.size());
  c.col_ind.reserve(a.values.size());

  std::ranges::fill(marker_, -1);
  used_cols_.clear();

  for (int i = 0; i < a.rows; ++i) {
    for (int j = a.row_ptr[i]; j < a.row_ptr[i + 1]; ++j) {
      const int col_a = a.col_ind[j];
      const double val_a = a.values[j];

      for (int k = b.row_ptr[col_a]; k < b.row_ptr[col_a + 1]; ++k) {
        const int col_b = b.col_ind[k];
        const double val_b = b.values[k];

        if (marker_[col_b] < i) {
          marker_[col_b] = i;
          used_cols_.push_back(col_b);
          temp_row_[col_b] = val_a * val_b;
        } else {
          temp_row_[col_b] += val_a * val_b;
        }
      }
    }

    if (!used_cols_.empty()) {
      std::ranges::sort(used_cols_);

      for (int col_idx : used_cols_) {
        const double val = temp_row_[col_idx];
        if (std::abs(val) > 1e-15) {
          c.values.push_back(val);
          c.col_ind.push_back(col_idx);
        }
        temp_row_[col_idx] = 0.0;
      }
      used_cols_.clear();
    }
    c.row_ptr[i + 1] = static_cast<int>(c.values.size());
  }

  return true;
}

bool MaslovaUMultMatrSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace maslova_u_mult_matr_crs
