#include "ermakov_a_spar_mat_mult_omp/omp/include/ops_omp.hpp"

#include <algorithm>
#include <complex>
#include <cstddef>
#include <vector>

#include "ermakov_a_spar_mat_mult_omp/common/include/common.hpp"

namespace ermakov_a_spar_mat_mult_omp {

ErmakovASparMatMultOMP::ErmakovASparMatMultOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool ErmakovASparMatMultOMP::ValidateMatrix(const MatrixCRS &m) {
  if (m.rows < 0 || m.cols < 0) {
    return false;
  }

  if (m.row_ptr.size() != static_cast<std::size_t>(m.rows) + 1) {
    return false;
  }
  if (m.values.size() != m.col_index.size()) {
    return false;
  }

  const int nnz = static_cast<int>(m.values.size());

  if (m.row_ptr.empty()) {
    return false;
  }
  if (m.row_ptr.front() != 0 || m.row_ptr.back() != nnz) {
    return false;
  }

  for (int i = 0; i < m.rows; ++i) {
    if (m.row_ptr[i] > m.row_ptr[i + 1]) {
      return false;
    }
  }

  for (int k = 0; k < nnz; ++k) {
    if (m.col_index[k] < 0 || m.col_index[k] >= m.cols) {
      return false;
    }
  }

  return true;
}

bool ErmakovASparMatMultOMP::ValidationImpl() {
  const auto &a = GetInput().A;
  const auto &b = GetInput().B;

  if (a.cols != b.rows) {
    return false;
  }
  if (!ValidateMatrix(a)) {
    return false;
  }
  if (!ValidateMatrix(b)) {
    return false;
  }

  return true;
}

bool ErmakovASparMatMultOMP::PreProcessingImpl() {
  a_ = GetInput().A;
  b_ = GetInput().B;

  c_.rows = a_.rows;
  c_.cols = b_.cols;

  c_.values.clear();
  c_.col_index.clear();
  c_.row_ptr.assign(static_cast<std::size_t>(c_.rows) + 1, 0);

  return true;
}

void ErmakovASparMatMultOMP::ProcessRow(int i, std::vector<std::complex<double>> &row_vals, std::vector<int> &row_mark,
                                        std::vector<int> &used_cols,
                                        std::vector<std::vector<std::complex<double>>> &row_values,
                                        std::vector<std::vector<int>> &row_cols) {
  used_cols.clear();

  const int a_start = a_.row_ptr[i];
  const int a_end = a_.row_ptr[i + 1];

  for (int ak = a_start; ak < a_end; ++ak) {
    const int j = a_.col_index[ak];
    const auto a_ij = a_.values[ak];

    const int b_start = b_.row_ptr[j];
    const int b_end = b_.row_ptr[j + 1];

    for (int bk = b_start; bk < b_end; ++bk) {
      const int k = b_.col_index[bk];
      const auto b_jk = b_.values[bk];

      if (row_mark[k] != i) {
        row_mark[k] = i;
        row_vals[k] = a_ij * b_jk;
        used_cols.push_back(k);
      } else {
        row_vals[k] += a_ij * b_jk;
      }
    }
  }

  std::ranges::sort(used_cols);

  auto &cols = row_cols[static_cast<std::size_t>(i)];
  auto &vals = row_values[static_cast<std::size_t>(i)];

  cols.reserve(used_cols.size());
  vals.reserve(used_cols.size());

  for (int k : used_cols) {
    const auto v = row_vals[k];
    if (v != std::complex<double>(0.0, 0.0)) {
      cols.push_back(k);
      vals.push_back(v);
    }
  }
}

bool ErmakovASparMatMultOMP::RunImpl() {
  const int m = a_.rows;
  const int p = b_.cols;

  if (a_.cols != b_.rows) {
    return false;
  }

  c_.values.clear();
  c_.col_index.clear();
  std::ranges::fill(c_.row_ptr, 0);

  std::vector<std::vector<std::complex<double>>> row_values(static_cast<std::size_t>(m));
  std::vector<std::vector<int>> row_cols(static_cast<std::size_t>(m));

#pragma omp parallel default(none) shared(row_values, row_cols, m, p)
  {
    std::vector<std::complex<double>> row_vals(static_cast<std::size_t>(p), std::complex<double>(0.0, 0.0));
    std::vector<int> row_mark(static_cast<std::size_t>(p), -1);
    std::vector<int> used_cols;
    used_cols.reserve(256);

#pragma omp for schedule(static)
    for (int i = 0; i < m; ++i) {
      ProcessRow(i, row_vals, row_mark, used_cols, row_values, row_cols);
    }
  }

  int nnz = 0;
  for (int i = 0; i < m; ++i) {
    c_.row_ptr[i] = nnz;
    nnz += static_cast<int>(row_values[static_cast<std::size_t>(i)].size());
  }
  c_.row_ptr[m] = nnz;

  c_.values.reserve(static_cast<std::size_t>(nnz));
  c_.col_index.reserve(static_cast<std::size_t>(nnz));

  for (int i = 0; i < m; ++i) {
    const auto idx = static_cast<std::size_t>(i);

    c_.values.insert(c_.values.end(), row_values[idx].begin(), row_values[idx].end());

    c_.col_index.insert(c_.col_index.end(), row_cols[idx].begin(), row_cols[idx].end());
  }

  return true;
}

bool ErmakovASparMatMultOMP::PostProcessingImpl() {
  GetOutput() = c_;
  return true;
}

}  // namespace ermakov_a_spar_mat_mult_omp
