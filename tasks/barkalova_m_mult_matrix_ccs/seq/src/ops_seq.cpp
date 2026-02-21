#include "barkalova_m_mult_matrix_ccs/seq/include/ops_seq.hpp"

#include <cmath>
#include <complex>
#include <cstddef>
#include <exception>
#include <vector>

#include "barkalova_m_mult_matrix_ccs/common/include/common.hpp"

namespace barkalova_m_mult_matrix_ccs {

namespace {
constexpr double kEpsilon = 1e-10;
}  // namespace

BarkalovaMMultMatrixCcsSEQ::BarkalovaMMultMatrixCcsSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = CCSMatrix{};
}

bool BarkalovaMMultMatrixCcsSEQ::ValidationImpl() {
  const auto &[A, B] = GetInput();

  // Самое важное: можно ли умножать матрицы
  if (A.cols != B.rows) {
    return false;
  }

  // Проверка, что матрицы не пустые
  if (A.rows <= 0 || A.cols <= 0 || B.rows <= 0 || B.cols <= 0) {
    return false;
  }

  // Проверка структуры CCS (самое минимальное)
  if (A.col_ptrs.size() != static_cast<size_t>(A.cols) + 1 || B.col_ptrs.size() != static_cast<size_t>(B.cols) + 1) {
    return false;
  }

  // Проверка, что col_ptrs начинается с 0
  if (A.col_ptrs.empty() || A.col_ptrs[0] != 0 || B.col_ptrs.empty() || B.col_ptrs[0] != 0) {
    return false;
  }

  // Проверка соответствия nnz
  if (A.nnz != static_cast<int>(A.values.size()) || B.nnz != static_cast<int>(B.values.size())) {
    return false;
  }

  return true;
}

bool BarkalovaMMultMatrixCcsSEQ::PreProcessingImpl() {
  return true;
}

void BarkalovaMMultMatrixCcsSEQ::TransposeMatrix(const CCSMatrix &a, CCSMatrix &at) {
  at.rows = a.cols;
  at.cols = a.rows;
  at.nnz = a.nnz;

  if (a.nnz == 0) {
    at.values.clear();
    at.row_indices.clear();
    at.col_ptrs.assign(at.cols + 1, 0);
    return;
  }

  std::vector<int> row_count(at.cols, 0);
  for (int i = 0; i < a.nnz; i++) {
    row_count[a.row_indices[i]]++;
  }

  at.col_ptrs.resize(at.cols + 1);
  at.col_ptrs[0] = 0;
  for (int i = 0; i < at.cols; i++) {
    at.col_ptrs[i + 1] = at.col_ptrs[i] + row_count[i];
  }

  at.values.resize(a.nnz);
  at.row_indices.resize(a.nnz);

  std::vector<int> current_pos(at.cols, 0);
  for (int col = 0; col < a.cols; col++) {
    for (int i = a.col_ptrs[col]; i < a.col_ptrs[col + 1]; i++) {
      int row = a.row_indices[i];
      Complex val = a.values[i];

      int pos = at.col_ptrs[row] + current_pos[row];
      at.values[pos] = val;
      at.row_indices[pos] = col;
      current_pos[row]++;
    }
  }
}

namespace {

void ProcessColumnSEQ(const CCSMatrix &at, const CCSMatrix &b, int col_index, std::vector<Complex> &temp_row,
                      std::vector<int> &row_marker, std::vector<Complex> &res_val, std::vector<int> &res_row_ind) {
  for (int k = b.col_ptrs[col_index]; k < b.col_ptrs[col_index + 1]; k++) {
    int row_b = b.row_indices[k];
    Complex val_b = b.values[k];

    for (int idx = at.col_ptrs[row_b]; idx < at.col_ptrs[row_b + 1]; idx++) {
      int row_a = at.row_indices[idx];
      Complex val_a = at.values[idx];

      if (row_marker[row_a] != col_index) {
        row_marker[row_a] = col_index;
        temp_row[row_a] = val_a * val_b;
      } else {
        temp_row[row_a] += val_a * val_b;
      }
    }
  }

  for (size_t i = 0; i < temp_row.size(); i++) {
    if (row_marker[i] == col_index && std::norm(temp_row[i]) > kEpsilon) {
      res_val.push_back(temp_row[i]);
      res_row_ind.push_back(static_cast<int>(i));
    }
  }
}

}  // namespace

void BarkalovaMMultMatrixCcsSEQ::MultiplyMatrices(const CCSMatrix &a, const CCSMatrix &b, CCSMatrix &c) {
  CCSMatrix at;
  TransposeMatrix(a, at);

  c.rows = a.rows;
  c.cols = b.cols;
  c.col_ptrs.push_back(0);

  std::vector<Complex> temp_row(c.rows, Complex(0.0, 0.0));
  std::vector<int> row_marker(c.rows, -1);

  for (int j = 0; j < b.cols; j++) {
    ProcessColumnSEQ(at, b, j, temp_row, row_marker, c.values, c.row_indices);
    c.col_ptrs.push_back(static_cast<int>(c.values.size()));
  }

  c.nnz = static_cast<int>(c.values.size());
}

bool BarkalovaMMultMatrixCcsSEQ::RunImpl() {
  const auto &[a, b] = GetInput();

  try {
    CCSMatrix c;
    MultiplyMatrices(a, b, c);
    GetOutput() = c;
    return true;
  } catch (const std::exception &) {
    return false;
  }
}

bool BarkalovaMMultMatrixCcsSEQ::PostProcessingImpl() {
  const auto &c = GetOutput();
  return c.rows > 0 && c.cols > 0 && c.col_ptrs.size() == static_cast<size_t>(c.cols) + 1;
}

}  // namespace barkalova_m_mult_matrix_ccs
