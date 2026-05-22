#include "borunov_v_complex_ccs/all/include/ops_all.hpp"

#include <omp.h>
#include <tbb/tbb.h>

#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>

#include "borunov_v_complex_ccs/common/include/common.hpp"
#include "util/include/util.hpp"

namespace borunov_v_complex_ccs {

BorunovVComplexCcsALL::BorunovVComplexCcsALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().resize(1);
}

bool BorunovVComplexCcsALL::ValidationImpl() {
  const auto &a = GetInput().first;
  const auto &b = GetInput().second;
  if (a.num_cols != b.num_rows) {
    return false;
  }
  if (a.col_ptrs.size() != static_cast<std::size_t>(a.num_cols) + 1 ||
      b.col_ptrs.size() != static_cast<std::size_t>(b.num_cols) + 1) {
    return false;
  }
  return true;
}

bool BorunovVComplexCcsALL::PreProcessingImpl() {
  const auto &a = GetInput().first;
  const auto &b = GetInput().second;
  auto &c = GetOutput()[0];

  c.num_rows = a.num_rows;
  c.num_cols = b.num_cols;
  c.col_ptrs.assign(c.num_cols + 1, 0);
  c.values.clear();
  c.row_indices.clear();

  return true;
}

namespace {

void CustomShellSort(std::vector<int> &touched) {
  for (std::size_t gap = touched.size() / 2; gap > 0; gap /= 2) {
    for (std::size_t i = gap; i < touched.size(); ++i) {
      const int temp = touched[i];
      std::size_t j = i;
      while (j >= gap && touched[j - gap] > temp) {
        touched[j] = touched[j - gap];
        j -= gap;
      }
      touched[j] = temp;
    }
  }
}

void ProcessColumn(int j, const SparseMatrix &a, const SparseMatrix &b, int local_j,
                   std::vector<std::complex<double>> &acc, std::vector<int> &marker, std::vector<int> &touched,
                   std::vector<std::complex<double>> &t_vals, std::vector<int> &t_rows, std::vector<int> &t_col_nnz) {
  touched.clear();

  for (int bk = b.col_ptrs[j]; bk < b.col_ptrs[j + 1]; ++bk) {
    const int p = b.row_indices[bk];
    const std::complex<double> bval = b.values[bk];

    for (int ak = a.col_ptrs[p]; ak < a.col_ptrs[p + 1]; ++ak) {
      const int i = a.row_indices[ak];
      acc[i] += a.values[ak] * bval;
      if (marker[i] != j) {
        marker[i] = j;
        touched.push_back(i);
      }
    }
  }

  CustomShellSort(touched);

  for (int i : touched) {
    if (std::abs(acc[i]) > 1e-9) {
      t_vals.push_back(acc[i]);
      t_rows.push_back(i);
      ++t_col_nnz[local_j];
    }
    acc[i] = {0.0, 0.0};
  }
}

void ComputeTaskRange(int start_col, int end_col, int tid, const SparseMatrix &a, const SparseMatrix &b,
                      std::vector<std::vector<std::complex<double>>> &vals, std::vector<std::vector<int>> &rows,
                      std::vector<std::vector<int>> &col_nnz) {
  const int cols_count = end_col - start_col;
  col_nnz[tid].assign(cols_count, 0);

  std::vector<std::complex<double>> acc(a.num_rows, {0.0, 0.0});
  std::vector<int> marker(a.num_rows, -1);
  std::vector<int> touched;
  touched.reserve(static_cast<std::size_t>(a.num_rows));

  for (int j = start_col; j < end_col; ++j) {
    ProcessColumn(j, a, b, j - start_col, acc, marker, touched, vals[tid], rows[tid], col_nnz[tid]);
  }
}

void FillColumnPointers(int start_col, int end_col, int num_threads, const std::vector<std::vector<int>> &col_nnz,
                        std::vector<int> &col_ptrs) {
  const int range = end_col - start_col;
  for (int tid = 0; tid < num_threads; ++tid) {
    const int tid_start_col = start_col + ((range * tid) / num_threads);
    const int tid_end_col = start_col + ((range * (tid + 1)) / num_threads);
    for (int j = tid_start_col; j < tid_end_col; ++j) {
      col_ptrs[j + 1] = col_ptrs[j] + col_nnz[tid][j - tid_start_col];
    }
  }
}

}  // namespace

bool BorunovVComplexCcsALL::RunImpl() {
  const auto &a = GetInput().first;
  const auto &b = GetInput().second;
  auto &c = GetOutput()[0];

  const int num_threads = ppc::util::GetNumThreads();
  const int num_cols = b.num_cols;

  if (num_cols == 0) {
    return true;
  }

  const int half_cols = num_cols / 2;

  std::vector<std::vector<std::complex<double>>> omp_vals(num_threads);
  std::vector<std::vector<int>> omp_rows(num_threads);
  std::vector<std::vector<int>> omp_col_nnz(num_threads);

  if (half_cols > 0) {
#pragma omp parallel num_threads(num_threads) default(none) \
    shared(a, b, half_cols, num_threads, omp_vals, omp_rows, omp_col_nnz)
    {
      const int tid = omp_get_thread_num();
      const int total_threads = omp_get_num_threads();
      const int tid_start_col = (half_cols * tid) / total_threads;
      const int tid_end_col = (half_cols * (tid + 1)) / total_threads;
      ComputeTaskRange(tid_start_col, tid_end_col, tid, a, b, omp_vals, omp_rows, omp_col_nnz);
    }
  }

  FillColumnPointers(0, half_cols, num_threads, omp_col_nnz, c.col_ptrs);

  std::vector<std::vector<std::complex<double>>> tbb_vals(num_threads);
  std::vector<std::vector<int>> tbb_rows(num_threads);
  std::vector<std::vector<int>> tbb_col_nnz(num_threads);
  if (num_cols > half_cols) {
    tbb::task_arena arena(num_threads);
    arena.execute([&] {
      tbb::parallel_for(tbb::blocked_range<int>(0, num_threads, 1), [&](const tbb::blocked_range<int> &r) {
        for (int tid = r.begin(); tid < r.end(); ++tid) {
          const int tid_start_col = half_cols + (((num_cols - half_cols) * tid) / num_threads);
          const int tid_end_col = half_cols + (((num_cols - half_cols) * (tid + 1)) / num_threads);
          ComputeTaskRange(tid_start_col, tid_end_col, tid, a, b, tbb_vals, tbb_rows, tbb_col_nnz);
        }
      }, tbb::static_partitioner());
    });
  }

  FillColumnPointers(half_cols, num_cols, num_threads, tbb_col_nnz, c.col_ptrs);

  int total_nnz = c.col_ptrs[num_cols];
  c.values.reserve(total_nnz);
  c.row_indices.reserve(total_nnz);

  if (half_cols > 0) {
    for (int tid = 0; tid < num_threads; ++tid) {
      c.values.insert(c.values.end(), omp_vals[tid].begin(), omp_vals[tid].end());
      c.row_indices.insert(c.row_indices.end(), omp_rows[tid].begin(), omp_rows[tid].end());
    }
  }

  if (num_cols > half_cols) {
    for (int tid = 0; tid < num_threads; ++tid) {
      c.values.insert(c.values.end(), tbb_vals[tid].begin(), tbb_vals[tid].end());
      c.row_indices.insert(c.row_indices.end(), tbb_rows[tid].begin(), tbb_rows[tid].end());
    }
  }

  return true;
}

bool BorunovVComplexCcsALL::PostProcessingImpl() {
  return true;
}

}  // namespace borunov_v_complex_ccs
