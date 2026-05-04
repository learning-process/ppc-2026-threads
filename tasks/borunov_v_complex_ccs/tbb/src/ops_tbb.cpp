#include "borunov_v_complex_ccs/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>

#include "borunov_v_complex_ccs/common/include/common.hpp"
#include "util/include/util.hpp"

namespace borunov_v_complex_ccs {

BorunovVComplexCcsTBB::BorunovVComplexCcsTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().resize(1);
}

bool BorunovVComplexCcsTBB::ValidationImpl() {
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

bool BorunovVComplexCcsTBB::PreProcessingImpl() {
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

void ProcessColumn(int j, const SparseMatrix &a, const SparseMatrix &b, int tid, int jstart,
                   std::vector<std::complex<double>> &acc, std::vector<int> &marker, std::vector<int> &touched,
                   std::vector<std::vector<std::complex<double>>> &t_values,
                   std::vector<std::vector<int>> &t_row_indices, std::vector<std::vector<int>> &t_col_nnz) {
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
      t_values[tid].push_back(acc[i]);
      t_row_indices[tid].push_back(i);
      ++t_col_nnz[tid][j - jstart];
    }
    acc[i] = {0.0, 0.0};
  }
}

}  // namespace

bool BorunovVComplexCcsTBB::RunImpl() {
  const auto &a = GetInput().first;
  const auto &b = GetInput().second;
  auto &c = GetOutput()[0];

  const int num_threads = ppc::util::GetNumThreads();
  const int bc = b.num_cols;

  std::vector<std::vector<std::complex<double>>> t_values(num_threads);
  std::vector<std::vector<int>> t_row_indices(num_threads);
  std::vector<std::vector<int>> t_col_nnz(num_threads);

  tbb::task_arena arena(num_threads);
  arena.execute([&] {
    tbb::parallel_for(tbb::blocked_range<int>(0, num_threads, 1), [&](const tbb::blocked_range<int> &r) {
      for (int tid = r.begin(); tid < r.end(); ++tid) {
        const int jstart = (tid * bc) / num_threads;
        const int jend = ((tid + 1) * bc) / num_threads;

        t_col_nnz[tid].assign(jend - jstart, 0);

        std::vector<std::complex<double>> acc(a.num_rows, {0.0, 0.0});
        std::vector<int> marker(a.num_rows, -1);
        std::vector<int> touched;
        touched.reserve(static_cast<std::size_t>(a.num_rows));

        for (int j = jstart; j < jend; ++j) {
          ProcessColumn(j, a, b, tid, jstart, acc, marker, touched, t_values, t_row_indices, t_col_nnz);
        }
      }
    }, tbb::static_partitioner());
  });

  for (int tid = 0; tid < num_threads; ++tid) {
    const int jstart = (tid * bc) / num_threads;
    const int jend = ((tid + 1) * bc) / num_threads;
    for (int j = jstart; j < jend; ++j) {
      c.col_ptrs[j + 1] = c.col_ptrs[j] + t_col_nnz[tid][j - jstart];
    }
  }

  const int total_nnz = c.col_ptrs[bc];
  c.values.resize(static_cast<std::size_t>(total_nnz));
  c.row_indices.resize(static_cast<std::size_t>(total_nnz));

  std::vector<int> thread_offsets(num_threads + 1, 0);
  for (int tid = 0; tid < num_threads; ++tid) {
    thread_offsets[tid + 1] = thread_offsets[tid] + static_cast<int>(t_values[tid].size());
  }

  arena.execute([&] {
    tbb::parallel_for(tbb::blocked_range<int>(0, num_threads, 1), [&](const tbb::blocked_range<int> &r) {
      for (int tid = r.begin(); tid < r.end(); ++tid) {
        std::copy(t_values[tid].begin(), t_values[tid].end(), c.values.begin() + thread_offsets[tid]);
        std::copy(t_row_indices[tid].begin(), t_row_indices[tid].end(), c.row_indices.begin() + thread_offsets[tid]);
      }
    }, tbb::static_partitioner());
  });

  return true;
}

bool BorunovVComplexCcsTBB::PostProcessingImpl() {
  return true;
}

}  // namespace borunov_v_complex_ccs
