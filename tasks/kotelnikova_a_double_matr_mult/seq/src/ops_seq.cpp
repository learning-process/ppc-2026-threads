#include "kotelnikova_a_double_matr_mult/seq/include/ops_seq.hpp"

#include <cmath>
#include <cstddef>
#include <vector>

#include "kotelnikova_a_double_matr_mult/common/include/common.hpp"

namespace kotelnikova_a_double_matr_mult {

KotelnikovaATaskSEQ::KotelnikovaATaskSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = SparseMatrixCCS();
}

bool KotelnikovaATaskSEQ::IsMatrixValid(const SparseMatrixCCS &matrix) {
  if (matrix.rows < 0 || matrix.cols < 0) {
    return false;
  }
  if (matrix.col_ptrs.size() != static_cast<size_t>(matrix.cols) + 1) {
    return false;
  }
  if (matrix.values.size() != matrix.row_indices.size()) {
    return false;
  }

  if (matrix.col_ptrs.empty() || matrix.col_ptrs[0] != 0) {
    return false;
  }

  const int total_elements = static_cast<int>(matrix.values.size());
  if (matrix.col_ptrs[matrix.cols] != total_elements) {
    return false;
  }

  for (size_t i = 0; i < matrix.col_ptrs.size() - 1; ++i) {
    if (matrix.col_ptrs[i] > matrix.col_ptrs[i + 1] || matrix.col_ptrs[i] < 0) {
      return false;
    }
  }

  for (size_t i = 0; i < matrix.row_indices.size(); ++i) {
    if (matrix.row_indices[i] < 0 || matrix.row_indices[i] >= matrix.rows) {
      return false;
    }
  }

  return true;
}

bool KotelnikovaATaskSEQ::ValidationImpl() {
  const auto &[a, b] = GetInput();

  if (!IsMatrixValid(a) || !IsMatrixValid(b)) {
    return false;
  }
  if (a.cols != b.rows) {
    return false;
  }

  return true;
}

bool KotelnikovaATaskSEQ::PreProcessingImpl() {
  const auto &[a, b] = GetInput();
  GetOutput() = SparseMatrixCCS(a.rows, b.cols);
  return true;
}

namespace {
void ProcessColumn(const SparseMatrixCCS &a, const SparseMatrixCCS &b, int j, std::vector<double> &temp) {
  for (int k = 0; k < a.cols; ++k) {
    double b_val = 0.0;
    for (int b_idx = b.col_ptrs[j]; b_idx < b.col_ptrs[j + 1]; ++b_idx) {
      if (b.row_indices[b_idx] == k) {
        b_val = b.values[b_idx];
        break;
      }
    }

    if (b_val == 0.0) {
      continue;
    }

    for (int a_idx = a.col_ptrs[k]; a_idx < a.col_ptrs[k + 1]; ++a_idx) {
      const int i = a.row_indices[a_idx];
      const double a_val = a.values[a_idx];
      temp[i] += a_val * b_val;
    }
  }
}

void SaveColumnResults(int rows, std::vector<double> &temp, SparseMatrixCCS &result) {
  for (int i = 0; i < rows; ++i) {
    if (std::abs(temp[i]) > 1e-10) {
      result.values.push_back(temp[i]);
      result.row_indices.push_back(i);
      temp[i] = 0.0;
    }
  }
}
}  // namespace

SparseMatrixCCS KotelnikovaATaskSEQ::MultiplyMatrices(const SparseMatrixCCS &a, const SparseMatrixCCS &b) {
  SparseMatrixCCS result(a.rows, b.cols);
  std::vector<double> temp(a.rows, 0.0);

  for (int j = 0; j < b.cols; ++j) {
    result.col_ptrs[j] = static_cast<int>(result.values.size());
    ProcessColumn(a, b, j, temp);
    SaveColumnResults(a.rows, temp, result);
  }

  result.col_ptrs[b.cols] = static_cast<int>(result.values.size());
  return result;
}

bool KotelnikovaATaskSEQ::RunImpl() {
  const auto &[a, b] = GetInput();
  GetOutput() = MultiplyMatrices(a, b);
  return true;
}

bool KotelnikovaATaskSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace kotelnikova_a_double_matr_mult
