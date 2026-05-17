#include "kulik_a_mat_mul_double_ccs/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cstddef>
#include <functional>
#include <thread>
#include <tuple>
#include <vector>

#include "kulik_a_mat_mul_double_ccs/common/include/common.hpp"

namespace kulik_a_mat_mul_double_ccs {

namespace {

inline void ProcessColumn(size_t j, const CCS &a, const CCS &b, std::vector<double> &accum,
                          std::vector<bool> &nz_elem_rows, std::vector<size_t> &nnz_rows,
                          std::vector<std::vector<double>> &local_values,
                          std::vector<std::vector<size_t>> &local_rows) {
  for (size_t k = b.col_ind[j]; k < b.col_ind[j + 1]; ++k) {
    size_t ind = b.row[k];
    double b_val = b.value[k];
    for (size_t zc = a.col_ind[ind]; zc < a.col_ind[ind + 1]; ++zc) {
      size_t i = a.row[zc];
      double a_val = a.value[zc];

      accum[i] += a_val * b_val;
      if (!nz_elem_rows[i]) {
        nz_elem_rows[i] = true;
        nnz_rows.push_back(i);
      }
    }
  }

  std::ranges::sort(nnz_rows);

  for (size_t i : nnz_rows) {
    if (accum[i] != 0.0) {
      local_rows[j].push_back(i);
      local_values[j].push_back(accum[i]);
    }
    accum[i] = 0.0;
    nz_elem_rows[i] = false;
  }
  nnz_rows.clear();
}

inline void CopyColumn(size_t j, CCS &c, const std::vector<std::vector<double>> &local_values,
                       const std::vector<std::vector<size_t>> &local_rows) {
  size_t offset = c.col_ind[j];
  size_t col_nz = local_values[j].size();
  for (size_t k = 0; k < col_nz; ++k) {
    c.value[offset + k] = local_values[j][k];
    c.row[offset + k] = local_rows[j][k];
  }
}

void ProcessColumnsRange(size_t jstart, size_t jend, const CCS &a, const CCS &b,
                         std::vector<std::vector<double>> &local_values, std::vector<std::vector<size_t>> &local_rows) {
  std::vector<double> accum(a.n, 0.0);
  std::vector<bool> nz_elem_rows(a.n, false);
  std::vector<size_t> nnz_rows;
  nnz_rows.reserve(a.n);

  for (size_t j = jstart; j < jend; ++j) {
    ProcessColumn(j, a, b, accum, nz_elem_rows, nnz_rows, local_values, local_rows);
  }
}

void CopyColumnsRange(size_t jstart, size_t jend, CCS &c, const std::vector<std::vector<double>> &local_values,
                      const std::vector<std::vector<size_t>> &local_rows) {
  for (size_t j = jstart; j < jend; ++j) {
    CopyColumn(j, c, local_values, local_rows);
  }
}

}  // namespace

KulikAMatMulDoubleCcsSTL::KulikAMatMulDoubleCcsSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool KulikAMatMulDoubleCcsSTL::ValidationImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());
  return (a.m == b.n);
}

bool KulikAMatMulDoubleCcsSTL::PreProcessingImpl() {
  return true;
}

bool KulikAMatMulDoubleCcsSTL::RunImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());
  OutType &c = GetOutput();

  c.n = a.n;
  c.m = b.m;
  c.col_ind.assign(c.m + 1, 0);

  std::vector<std::vector<double>> local_values(b.m);
  std::vector<std::vector<size_t>> local_rows(b.m);

  const size_t num_threads_raw = std::thread::hardware_concurrency();
  const size_t num_threads = std::max<size_t>(1, num_threads_raw == 0 ? 1 : num_threads_raw);
  const size_t threads_count = std::min(num_threads, b.m == 0 ? static_cast<size_t>(1) : b.m);

  std::vector<std::thread> threads;
  threads.reserve(threads_count);

  for (size_t tid = 0; tid < threads_count; ++tid) {
    const size_t jstart = (tid * b.m) / threads_count;
    const size_t jend = ((tid + 1) * b.m) / threads_count;

    threads.emplace_back(ProcessColumnsRange, jstart, jend, std::cref(a), std::cref(b), std::ref(local_values),
                         std::ref(local_rows));
  }

  for (auto &th : threads) {
    th.join();
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

  threads.clear();
  for (size_t tid = 0; tid < threads_count; ++tid) {
    const size_t jstart = (tid * b.m) / threads_count;
    const size_t jend = ((tid + 1) * b.m) / threads_count;

    threads.emplace_back(CopyColumnsRange, jstart, jend, std::ref(c), std::cref(local_values), std::cref(local_rows));
  }

  for (auto &th : threads) {
    th.join();
  }

  return true;
}

bool KulikAMatMulDoubleCcsSTL::PostProcessingImpl() {
  return true;
}

}  // namespace kulik_a_mat_mul_double_ccs
