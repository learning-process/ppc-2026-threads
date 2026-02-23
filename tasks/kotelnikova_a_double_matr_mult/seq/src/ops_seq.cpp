#include "kotelnikova_a_double_matr_mult/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace kotelnikova_a_double_matr_mult {

KotelnikovaATaskSEQ::KotelnikovaATaskSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = SparseMatrixCCS();
}

bool KotelnikovaATaskSEQ::IsMatrixValid(const SparseMatrixCCS &matrix) {
  if (matrix.rows < 0 || matrix.cols < 0) return false;
  if (matrix.col_ptrs.size() != static_cast<size_t>(matrix.cols + 1)) return false;
  if (matrix.values.size() != matrix.row_indices.size()) return false;
  
  for (size_t i = 0; i < matrix.col_ptrs.size() - 1; ++i) {
    if (matrix.col_ptrs[i] > matrix.col_ptrs[i + 1]) return false;
    if (matrix.col_ptrs[i] < 0 || matrix.col_ptrs[i + 1] < 0) return false;
  }
  
  int total_elements = static_cast<int>(matrix.values.size());
  if (matrix.col_ptrs[0] != 0) return false;
  if (matrix.col_ptrs[matrix.cols] != total_elements) return false;
  
  for (size_t i = 0; i < matrix.row_indices.size(); ++i) {
    if (matrix.row_indices[i] < 0 || matrix.row_indices[i] >= matrix.rows) return false;
  }
  
  return true;
}

bool KotelnikovaATaskSEQ::ValidationImpl() {
  const auto &[A, B] = GetInput();
  
  if (!IsMatrixValid(A) || !IsMatrixValid(B)) return false;
  if (A.cols != B.rows) return false;
  
  return true;
}

bool KotelnikovaATaskSEQ::PreProcessingImpl() {
  const auto &[A, B] = GetInput();
  GetOutput() = SparseMatrixCCS(A.rows, B.cols);
  return true;
}

SparseMatrixCCS KotelnikovaATaskSEQ::MultiplyMatrices(const SparseMatrixCCS &A, const SparseMatrixCCS &B) {
  SparseMatrixCCS result(A.rows, B.cols);
  
  std::vector<double> temp(A.rows, 0.0);
  
  for (int j = 0; j < B.cols; ++j) {
    result.col_ptrs[j] = static_cast<int>(result.values.size());
    
    for (int k = 0; k < A.cols; ++k) {
      double b_val = 0.0;
      for (int b_idx = B.col_ptrs[j]; b_idx < B.col_ptrs[j + 1]; ++b_idx) {
        if (B.row_indices[b_idx] == k) {
          b_val = B.values[b_idx];
          break;
        }
      }
      
      if (b_val == 0.0) continue;
      
      for (int a_idx = A.col_ptrs[k]; a_idx < A.col_ptrs[k + 1]; ++a_idx) {
        int i = A.row_indices[a_idx];
        double a_val = A.values[a_idx];
        temp[i] += a_val * b_val;
      }
    }
    
    for (int i = 0; i < A.rows; ++i) {
      if (std::abs(temp[i]) > 1e-10) {
        result.values.push_back(temp[i]);
        result.row_indices.push_back(i);
        temp[i] = 0.0;
      }
    }
  }
  
  result.col_ptrs[B.cols] = static_cast<int>(result.values.size());
  return result;
}

bool KotelnikovaATaskSEQ::RunImpl() {
  const auto &[A, B] = GetInput();
  GetOutput() = MultiplyMatrices(A, B);
  return true;
}

bool KotelnikovaATaskSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace kotelnikova_a_double_matr_mult