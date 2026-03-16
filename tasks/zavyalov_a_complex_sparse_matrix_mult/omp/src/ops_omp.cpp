#include "zavyalov_a_complex_sparse_matrix_mult/common/include/common.hpp"
#include "zavyalov_a_complex_sparse_matrix_mult/omp/include/ops_omp.hpp"

#include <atomic>
#include <numeric>
#include <vector>

#include "util/include/util.hpp"

namespace zavyalov_a_compl_sparse_matr_mult {

SparseMatrix ZavyalovAComplSparseMatrMultOMP::multiplicate_with_omp(const SparseMatrix &matr_a, const SparseMatrix &matr_b) {
  if (matr_a.width != matr_b.height) {
    throw std::invalid_argument("Incompatible matrix dimensions for multiplication");
  }

  std::map<std::pair<size_t, size_t>, Complex> mp;  // <row, col> -> val

  for (size_t i = 0; i < matr_a.Count(); ++i) {
    size_t row_a = matr_a.row_ind[i];
    size_t col_a = matr_a.col_ind[i];
    Complex val_a = matr_a.val[i];

    for (size_t j = 0; j < matr_b.Count(); ++j) {
      size_t row_b = matr_b.row_ind[j];
      size_t col_b = matr_b.col_ind[j];
      Complex val_b = matr_b.val[j];

      if (col_a == row_b) {
        mp[{row_a, col_b}] += val_a * val_b;
      }
    }
  }

  SparseMatrix res;
  res.width = matr_b.width;
  res.height = matr_a.height;
  for (const auto &[key, value] : mp) {
    res.val.push_back(value);
    res.row_ind.push_back(key.first);
    res.col_ind.push_back(key.second);
  }

  return res;
}

ZavyalovAComplSparseMatrMultOMP::ZavyalovAComplSparseMatrMultOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool ZavyalovAComplSparseMatrMultOMP::ValidationImpl() {
  const auto &matr_a = std::get<0>(GetInput());
  const auto &matr_b = std::get<1>(GetInput());
  return matr_a.width == matr_b.height;
}

bool ZavyalovAComplSparseMatrMultOMP::PreProcessingImpl() {
  return true;
}

bool ZavyalovAComplSparseMatrMultOMP::RunImpl() {
  const auto &matr_a = std::get<0>(GetInput());
  const auto &matr_b = std::get<1>(GetInput());


  GetOutput() = multiplicate_with_omp(matr_a, matr_b);

  return true;
}

bool ZavyalovAComplSparseMatrMultOMP::PostProcessingImpl() {
  return true;
}
}  // namespace zavyalov_a_compl_sparse_matr_mult
