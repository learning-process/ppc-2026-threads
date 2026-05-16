#include "maslova_u_mult_matr_crs/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "maslova_u_mult_matr_crs/common/include/common.hpp"
#include "util/include/util.hpp"

namespace maslova_u_mult_matr_crs {

MaslovaUMultMatrOMP::MaslovaUMultMatrOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool MaslovaUMultMatrOMP::ValidationImpl() {
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

bool MaslovaUMultMatrOMP::PreProcessingImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());
  auto &c = GetOutput();
  c.rows = a.rows;
  c.cols = b.cols;
  return true;
}

int MaslovaUMultMatrOMP::GetRowNNZ(int i, const CRSMatrix &a, const CRSMatrix &b, std::vector<int> &marker) {
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

void MaslovaUMultMatrOMP::FillRowValues(int i, const CRSMatrix &a, const CRSMatrix &b, CRSMatrix &c,
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

bool MaslovaUMultMatrOMP::RunImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());
  auto &c = GetOutput();

  const int rows_a = a.rows;
  const int cols_b = b.cols;
  c.row_ptr.assign(static_cast<size_t>(rows_a) + 1, 0);

#pragma omp parallel default(none) shared(a, b, c, rows_a, cols_b) num_threads(ppc::util::GetNumThreads())
  {
    std::vector<int> marker(cols_b, -1);
#pragma omp for schedule(dynamic, 20)
    for (int i = 0; i < rows_a; ++i) {
      c.row_ptr[i + 1] = GetRowNNZ(i, a, b, marker);
    }
  }

  for (int i = 0; i < rows_a; ++i) {
    c.row_ptr[i + 1] += c.row_ptr[i];
  }

  c.values.resize(c.row_ptr[rows_a]);
  c.col_ind.resize(c.row_ptr[rows_a]);

#pragma omp parallel default(none) shared(a, b, c, rows_a, cols_b) num_threads(ppc::util::GetNumThreads())
  {
    std::vector<double> acc(cols_b, 0.0);
    std::vector<int> marker(cols_b, -1);
    std::vector<int> used;
    used.reserve(cols_b);

#pragma omp for schedule(dynamic, 20)
    for (int i = 0; i < rows_a; ++i) {
      FillRowValues(i, a, b, c, acc, marker, used);
    }
  }

  return true;
}

bool MaslovaUMultMatrOMP::PostProcessingImpl() {
  return true;
}

}  // namespace maslova_u_mult_matr_crs
