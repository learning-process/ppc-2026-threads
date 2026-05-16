#include "zavyalov_a_complex_sparse_matrix_mult/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstddef>
#include <functional>
#include <map>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

#include "util/include/util.hpp"
#include "zavyalov_a_complex_sparse_matrix_mult/common/include/common.hpp"

namespace zavyalov_a_compl_sparse_matr_mult {  // comm for ci1234

void ZavyalovAComplSparseMatrMultSTL::Worker(
    int tid, std::size_t start, std::size_t end, const SparseMatrix &matr_a, const SparseMatrix &matr_b,
    std::vector<std::map<std::pair<std::size_t, std::size_t>, Complex>> &local_maps) {
  for (std::size_t i = start; i < end; ++i) {
    std::size_t row_a = matr_a.row_ind[i];
    std::size_t col_a = matr_a.col_ind[i];
    Complex val_a = matr_a.val[i];

    for (std::size_t j = 0; j < matr_b.Count(); ++j) {
      if (col_a == matr_b.row_ind[j]) {
        local_maps[tid][{row_a, matr_b.col_ind[j]}] += val_a * matr_b.val[j];
      }
    }
  }
}

SparseMatrix ZavyalovAComplSparseMatrMultSTL::MultiplicateWithStl(const SparseMatrix &matr_a,
                                                                  const SparseMatrix &matr_b) {
  if (matr_a.width != matr_b.height) {
    throw std::invalid_argument("Incompatible matrix dimensions for multiplication");
  }

  int num_threads = ppc::util::GetNumThreads();
  std::size_t total = matr_a.Count();

  std::vector<std::map<std::pair<std::size_t, std::size_t>, Complex>> local_maps(num_threads);

  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  std::size_t chunk = (total + num_threads - 1) / num_threads;
  for (int ti = 0; ti < num_threads; ++ti) {
    std::size_t start = ti * chunk;
    std::size_t end = std::min(start + chunk, total);
    if (start < total) {
      threads.emplace_back(Worker, ti, start, end, std::cref(matr_a), std::cref(matr_b), std::ref(local_maps));
    }
  }

  for (auto &th : threads) {
    th.join();
  }

  std::map<std::pair<std::size_t, std::size_t>, Complex> mp;
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

ZavyalovAComplSparseMatrMultSTL::ZavyalovAComplSparseMatrMultSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool ZavyalovAComplSparseMatrMultSTL::ValidationImpl() {
  const auto &matr_a = std::get<0>(GetInput());
  const auto &matr_b = std::get<1>(GetInput());
  return matr_a.width == matr_b.height;
}

bool ZavyalovAComplSparseMatrMultSTL::PreProcessingImpl() {
  return true;
}

bool ZavyalovAComplSparseMatrMultSTL::RunImpl() {
  const auto &matr_a = std::get<0>(GetInput());
  const auto &matr_b = std::get<1>(GetInput());

  GetOutput() = MultiplicateWithStl(matr_a, matr_b);

  return true;
}

bool ZavyalovAComplSparseMatrMultSTL::PostProcessingImpl() {
  return true;
}

}  // namespace zavyalov_a_compl_sparse_matr_mult
