#include "liulin_y_complex_ccs/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <atomic>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iterator>
#include <utility>
#include <vector>

#include "liulin_y_complex_ccs/common/include/common.hpp"

namespace liulin_y_complex_ccs {

namespace {

constexpr double kEpsilon = 1e-10;

bool IsValidCCS(const CCSMatrix &mat) {
  if (mat.count_rows <= 0 || mat.count_cols <= 0) {
    return false;
  }
  if (mat.col_index.size() != static_cast<size_t>(mat.count_cols) + 1) {
    return false;
  }
  if (mat.col_index[0] != 0) {
    return false;
  }
  if (mat.values.size() != mat.row_index.size()) {
    return false;
  }
  if (std::cmp_not_equal(mat.col_index.back(), mat.values.size())) {
    return false;
  }
  return true;
}

void TransposeCCS(const CCSMatrix &mat, CCSMatrix &mat_t) {
  mat_t.count_rows = mat.count_cols;
  mat_t.count_cols = mat.count_rows;
  const size_t nnz = mat.values.size();
  mat_t.values.resize(nnz);
  mat_t.row_index.resize(nnz);
  mat_t.col_index.assign(static_cast<size_t>(mat_t.count_cols) + 1, 0);

  for (size_t i = 0; i < nnz; ++i) {
    mat_t.col_index[static_cast<size_t>(mat.row_index[i]) + 1]++;
  }
  for (int i = 0; i < mat_t.count_cols; ++i) {
    mat_t.col_index[static_cast<size_t>(i) + 1] += mat_t.col_index[static_cast<size_t>(i)];
  }

  std::vector<int> current_pos(mat_t.col_index.begin(), mat_t.col_index.end());
  for (int j = 0; j < mat.count_cols; ++j) {
    for (int k = mat.col_index[static_cast<size_t>(j)]; k < mat.col_index[static_cast<size_t>(j) + 1]; ++k) {
      const int row = mat.row_index[static_cast<size_t>(k)];
      const int dest = current_pos[static_cast<size_t>(row)]++;
      mat_t.row_index[static_cast<size_t>(dest)] = j;
      mat_t.values[static_cast<size_t>(dest)] = mat.values[static_cast<size_t>(k)];
    }
  }
}

std::complex<double> ComputeDotProduct(int a_start, int a_end, int b_start, int b_end, const CCSMatrix &mat_at,
                                       const CCSMatrix &mat_b) {
  int ptr_a = a_start;
  int ptr_b = b_start;
  std::complex<double> sum(0.0, 0.0);

  while (ptr_a < a_end && ptr_b < b_end) {
    const int idx_a = mat_at.row_index[static_cast<size_t>(ptr_a)];
    const int idx_b = mat_b.row_index[static_cast<size_t>(ptr_b)];

    if (idx_a < idx_b) {
      ptr_a++;
    } else if (idx_a > idx_b) {
      ptr_b++;
    } else {
      sum += mat_at.values[static_cast<size_t>(ptr_a)] * mat_b.values[static_cast<size_t>(ptr_b)];
      ptr_a++;
      ptr_b++;
    }
  }
  return sum;
}

void ProcessColumn(int j, int res_rows, const CCSMatrix &mat_at, const CCSMatrix &mat_b,
                   std::vector<std::complex<double>> &col_values, std::vector<int> &col_rows) {
  const int b_start = mat_b.col_index[static_cast<size_t>(j)];
  const int b_end = mat_b.col_index[static_cast<size_t>(j) + 1];

  if (b_start == b_end) {
    return;
  }

  for (int i = 0; i < res_rows; ++i) {
    const int a_start = mat_at.col_index[static_cast<size_t>(i)];
    const int a_end = mat_at.col_index[static_cast<size_t>(i) + 1];

    if (a_start == a_end) {
      continue;
    }

    std::complex<double> sum = ComputeDotProduct(a_start, a_end, b_start, b_end, mat_at, mat_b);

    if (std::abs(sum.real()) > kEpsilon || std::abs(sum.imag()) > kEpsilon) {
      col_values.push_back(sum);
      col_rows.push_back(i);
    }
  }
}

}  // namespace

LiulinYComplexCcsTbb::LiulinYComplexCcsTbb(const InType &in) : BaseTask() {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool LiulinYComplexCcsTbb::ValidationImpl() {
  const auto &mat_a = GetInput().first;
  const auto &mat_b = GetInput().second;

  if (!IsValidCCS(mat_a) || !IsValidCCS(mat_b)) {
    return false;
  }
  if (mat_a.count_cols != mat_b.count_rows) {
    return false;
  }

  return true;
}

bool LiulinYComplexCcsTbb::PreProcessingImpl() {
  return true;
}

bool LiulinYComplexCcsTbb::RunImpl() {
  const auto &mat_a = GetInput().first;
  const auto &mat_b = GetInput().second;
  auto &mat_res = GetOutput();

  const int res_rows = mat_a.count_rows;
  const int res_cols = mat_b.count_cols;

  CCSMatrix mat_at;
  TransposeCCS(mat_a, mat_at);

  std::vector<std::vector<std::complex<double>>> thread_values(static_cast<size_t>(res_cols));
  std::vector<std::vector<int>> thread_row_indices(static_cast<size_t>(res_cols));
  std::atomic<bool> success{true};

  try {
    tbb::parallel_for(tbb::blocked_range<int>(0, res_cols), [&](const tbb::blocked_range<int> &r) {
      for (int j = r.begin(); j != r.end(); ++j) {
        ProcessColumn(j, res_rows, mat_at, mat_b, thread_values[static_cast<size_t>(j)],
                      thread_row_indices[static_cast<size_t>(j)]);
      }
    });
  } catch (...) {
    success = false;
  }

  if (!success) {
    return false;
  }

  mat_res.count_rows = res_rows;
  mat_res.count_cols = res_cols;
  mat_res.values.clear();
  mat_res.row_index.clear();
  mat_res.col_index.assign(static_cast<size_t>(res_cols) + 1, 0);

  for (int j = 0; j < res_cols; ++j) {
    const auto u_j = static_cast<size_t>(j);
    mat_res.values.insert(mat_res.values.end(), std::make_move_iterator(thread_values[u_j].begin()),
                          std::make_move_iterator(thread_values[u_j].end()));
    mat_res.row_index.insert(mat_res.row_index.end(), std::make_move_iterator(thread_row_indices[u_j].begin()),
                             std::make_move_iterator(thread_row_indices[u_j].end()));
    mat_res.col_index[u_j + 1] = static_cast<int>(mat_res.values.size());
  }

  return true;
}

bool LiulinYComplexCcsTbb::PostProcessingImpl() {
  const auto &mat_res = GetOutput();
  if (mat_res.count_rows <= 0 || mat_res.count_cols <= 0) {
    return false;
  }
  if (mat_res.col_index.size() != static_cast<size_t>(mat_res.count_cols) + 1) {
    return false;
  }
  return true;
}

}  // namespace liulin_y_complex_ccs
