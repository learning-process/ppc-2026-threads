#include "zavyalov_a_complex_sparse_matrix_mult/omp/include/ops_omp.hpp"

#include <omp.h>

#include <cstddef>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

#include "util/include/util.hpp"
#include "zavyalov_a_complex_sparse_matrix_mult/common/include/common.hpp"

namespace zavyalov_a_compl_sparse_matr_mult {

SparseMatrix ZavyalovAComplSparseMatrMultOMP::MultiplicateWithOmp(const SparseMatrix &matr_a,
                                                                  const SparseMatrix &matr_b) {
  if (matr_a.width != matr_b.height) {
    throw std::invalid_argument("Incompatible matrix dimensions for multiplication");
  }

  int num_threads = ppc::util::GetNumThreads();

  std::vector<std::map<std::pair<size_t, size_t>, Complex>> local_maps(num_threads);

#pragma omp parallel for num_threads(num_threads) schedule(static) default(none) shared(matr_a, matr_b, local_maps)
  for (size_t i = 0; i < matr_a.Count(); ++i) {
    int tid = omp_get_thread_num();
    size_t row_a = matr_a.row_ind[i];
    size_t col_a = matr_a.col_ind[i];
    Complex val_a = matr_a.val[i];

    for (size_t j = 0; j < matr_b.Count(); ++j) {
      if (col_a == matr_b.row_ind[j]) {
        local_maps[tid][{row_a, matr_b.col_ind[j]}] += val_a * matr_b.val[j];
      }
    }
  }

  std::map<std::pair<size_t, size_t>, Complex> mp;
  for (auto &lm : local_maps) {
    for (auto &[key, value] : lm) {
      mp[key] += value;
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

  GetOutput() = MultiplicateWithOmp(matr_a, matr_b);

  return true;
}

bool ZavyalovAComplSparseMatrMultOMP::PostProcessingImpl() {
  return true;
}
}  // namespace zavyalov_a_compl_sparse_matr_mult
