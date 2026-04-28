#include "zavyalov_a_complex_sparse_matrix_mult/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <cstddef>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

#include "zavyalov_a_complex_sparse_matrix_mult/common/include/common.hpp"

namespace zavyalov_a_compl_sparse_matr_mult {

SparseMatrix ZavyalovAComplSparseMatrMultTBB::MultiplicateWithTbb(const SparseMatrix &matr_a,
                                                                  const SparseMatrix &matr_b) {
  if (matr_a.width != matr_b.height) {
    throw std::invalid_argument("Incompatible matrix dimensions for multiplication");
  }

  tbb::enumerable_thread_specific<std::map<std::pair<size_t, size_t>, Complex>> local_maps;

  tbb::parallel_for(tbb::blocked_range<size_t>(0, matr_a.Count()), [&](const tbb::blocked_range<size_t> &range) {
    auto &my_map = local_maps.local();
    for (size_t i = range.begin(); i != range.end(); ++i) {
      size_t row_a = matr_a.row_ind[i];
      size_t col_a = matr_a.col_ind[i];
      Complex val_a = matr_a.val[i];

      for (size_t j = 0; j < matr_b.Count(); ++j) {
        if (col_a == matr_b.row_ind[j]) {
          my_map[{row_a, matr_b.col_ind[j]}] += val_a * matr_b.val[j];
        }
      }
    }
  });

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

ZavyalovAComplSparseMatrMultTBB::ZavyalovAComplSparseMatrMultTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool ZavyalovAComplSparseMatrMultTBB::ValidationImpl() {
  const auto &matr_a = std::get<0>(GetInput());
  const auto &matr_b = std::get<1>(GetInput());
  return matr_a.width == matr_b.height;
}

bool ZavyalovAComplSparseMatrMultTBB::PreProcessingImpl() {
  return true;
}

bool ZavyalovAComplSparseMatrMultTBB::RunImpl() {
  const auto &matr_a = std::get<0>(GetInput());
  const auto &matr_b = std::get<1>(GetInput());

  GetOutput() = MultiplicateWithTbb(matr_a, matr_b);

  return true;
}

bool ZavyalovAComplSparseMatrMultTBB::PostProcessingImpl() {
  return true;
}

}  // namespace zavyalov_a_compl_sparse_matr_mult
