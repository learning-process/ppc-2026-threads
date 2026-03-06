#include "kulik_a_mat_mul_double_ccs/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <tuple>
#include <vector>

#include "kulik_a_mat_mul_double_ccs/common/include/common.hpp"

namespace kulik_a_mat_mul_double_ccs {

void KulikAMatMulDoubleCcsOMP::ProcessColumn(size_t j, int tid, const CCS &a, const CCS &b,
                                             std::vector<std::vector<double>> &thread_accum,
                                             std::vector<std::vector<bool>> &thread_nz,
                                             std::vector<std::vector<size_t>> &thread_nnz_rows,
                                             std::vector<std::vector<double>> &local_values,
                                             std::vector<std::vector<size_t>> &local_rows) {
  for (size_t k = b.col_ind[j]; k < b.col_ind[j + 1]; ++k) {
    size_t ind = b.row[k];
    double b_val = b.value[k];
    for (size_t zc = a.col_ind[ind]; zc < a.col_ind[ind + 1]; ++zc) {
      size_t i = a.row[zc];
      double a_val = a.value[zc];
      thread_accum[tid][i] += a_val * b_val;
      if (!thread_nz[tid][i]) {
        thread_nz[tid][i] = true;
        thread_nnz_rows[tid].push_back(i);
      }
    }
  }

  std::ranges::sort(thread_nnz_rows[tid]);

  for (size_t i : thread_nnz_rows[tid]) {
    if (thread_accum[tid][i] != 0.0) {
      local_rows[j].push_back(i);
      local_values[j].push_back(thread_accum[tid][i]);
    }
    thread_accum[tid][i] = 0.0;
    thread_nz[tid][i] = false;
  }
  thread_nnz_rows[tid].clear();
}

KulikAMatMulDoubleCcsOMP::KulikAMatMulDoubleCcsOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool KulikAMatMulDoubleCcsOMP::ValidationImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());
  return (a.m == b.n);
}

bool KulikAMatMulDoubleCcsOMP::PreProcessingImpl() {
  return true;
}

bool KulikAMatMulDoubleCcsOMP::RunImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());
  OutType &c = GetOutput();
  c.n = a.n;
  c.m = b.m;
  c.col_ind.assign(c.m + 1, 0);

  std::vector<std::vector<double>> local_values(b.m);
  std::vector<std::vector<size_t>> local_rows(b.m);

  int num_threads = omp_get_max_threads();

  std::vector<std::vector<double>> thread_accum(num_threads, std::vector<double>(a.n, 0.0));
  std::vector<std::vector<bool>> thread_nz(num_threads, std::vector<bool>(a.n, false));
  std::vector<std::vector<size_t>> thread_nnz_rows(num_threads);

#pragma omp parallel for default(none) schedule(static) \
    shared(a, b, thread_accum, thread_nz, thread_nnz_rows, local_values, local_rows)
  for (size_t j = 0; j < b.m; ++j) {
    int tid = omp_get_thread_num();
    ProcessColumn(j, tid, a, b, thread_accum, thread_nz, thread_nnz_rows, local_values, local_rows);
  }

  size_t total_nz = 0;
  for (size_t j = 0; j < b.m; ++j) {
    c.col_ind[j] = total_nz;
    total_nz += local_values[j].size();
  }
  c.col_ind[b.m] = total_nz;
  c.nz = total_nz;

  c.value.resize(total_nz);
  c.row.resize(total_nz);

#pragma omp parallel for default(none) schedule(static) shared(b, c, local_values, local_rows)
  for (size_t j = 0; j < b.m; ++j) {
    size_t offset = c.col_ind[j];
    size_t col_nz = local_values[j].size();
    for (size_t k = 0; k < col_nz; ++k) {
      c.value[offset + k] = local_values[j][k];
      c.row[offset + k] = local_rows[j][k];
    }
  }

  return true;
}

bool KulikAMatMulDoubleCcsOMP::PostProcessingImpl() {
  return true;
}

}  // namespace kulik_a_mat_mul_double_ccs
