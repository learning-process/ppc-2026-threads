#include "maslova_u_mult_matr_crs/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <algorithm>
#include <vector>

#include "maslova_u_mult_matr_crs/common/include/common.hpp"

namespace maslova_u_mult_matr_crs {

MaslovaUMultMatrTBB::MaslovaUMultMatrTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool MaslovaUMultMatrTBB::ValidationImpl() {
  const auto &input = GetInput();
  const auto &a = std::get<0>(input);
  const auto &b = std::get<1>(input);
  return (a.cols == b.rows && a.rows > 0 && b.cols > 0);
}

bool MaslovaUMultMatrTBB::PreProcessingImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());
  auto &c = GetOutput();
  c.rows = a.rows;
  c.cols = b.cols;
  return true;
}

int MaslovaUMultMatrTBB::GetRowNNZ(int i, const CRSMatrix &a, const CRSMatrix &b, std::vector<int> &marker) {
  int row_nnz = 0;
  for (int j = a.row_ptr[i]; j < a.row_ptr[i + 1]; ++j) {
    int col_a = a.col_ind[j];
    for (int k = b.row_ptr[col_a]; k < b.row_ptr[col_a + 1]; ++k) {
      int col_b = b.col_ind[k];
      if (marker[col_b] != i) {
        marker[col_b] = i;
        row_nnz++;
      }
    }
  }
  return row_nnz;
}

void MaslovaUMultMatrTBB::FillRowValues(int i, const CRSMatrix &a, const CRSMatrix &b, CRSMatrix &c,
                                        std::vector<double> &acc, std::vector<int> &marker, std::vector<int> &used) {
  used.clear();
  for (int j = a.row_ptr[i]; j < a.row_ptr[i + 1]; ++j) {
    int col_a = a.col_ind[j];
    double val_a = a.values[j];
    for (int k = b.row_ptr[col_a]; k < b.row_ptr[col_a + 1]; ++k) {
      int col_b = b.col_ind[k];
      if (marker[col_b] != i) {
        marker[col_b] = i;
        used.push_back(col_b);
        acc[col_b] = val_a * b.values[k];
      } else {
        acc[col_b] += val_a * b.values[k];
      }
    }
  }
  std::ranges::sort(used);
  int write_pos = c.row_ptr[i];
  for (int col : used) {
    c.values[write_pos] = acc[col];
    c.col_ind[write_pos] = col;
    write_pos++;
    acc[col] = 0.0;
  }
}

bool MaslovaUMultMatrTBB::RunImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());
  auto &c = GetOutput();

  c.row_ptr.assign(a.rows + 1, 0);

  tbb::parallel_for(tbb::blocked_range<int>(0, a.rows), [&](const tbb::blocked_range<int> &r) {
    std::vector<int> marker(b.cols, -1);
    for (int i = r.begin(); i < r.end(); ++i) {
      c.row_ptr[i + 1] = GetRowNNZ(i, a, b, marker);
    }
  });

  for (int i = 0; i < a.rows; ++i) {
    c.row_ptr[i + 1] += c.row_ptr[i];
  }

  c.values.resize(c.row_ptr[a.rows]);
  c.col_ind.resize(c.row_ptr[a.rows]);

  tbb::parallel_for(tbb::blocked_range<int>(0, a.rows), [&](const tbb::blocked_range<int> &r) {
    std::vector<double> acc(b.cols, 0.0);
    std::vector<int> marker(b.cols, -1);
    std::vector<int> used;
    used.reserve(b.cols);
    for (int i = r.begin(); i < r.end(); ++i) {
      FillRowValues(i, a, b, c, acc, marker, used);
    }
  });

  return true;
}

bool MaslovaUMultMatrTBB::PostProcessingImpl() {
  return true;
}

}  // namespace maslova_u_mult_matr_crs
