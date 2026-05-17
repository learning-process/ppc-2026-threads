#include "barkalova_m_mult_matrix_ccs/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <exception>
#include <functional>
#include <future>
#include <thread>
#include <utility>
#include <vector>

#include "barkalova_m_mult_matrix_ccs/common/include/common.hpp"

namespace barkalova_m_mult_matrix_ccs {

BarkalovaMMultMatrixCcsSTL::BarkalovaMMultMatrixCcsSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = CCSMatrix{};
}

bool BarkalovaMMultMatrixCcsSTL::ValidationImpl() {
  const auto &[A, B] = GetInput();
  if (A.cols != B.rows) {
    return false;
  }
  if (A.rows <= 0 || A.cols <= 0 || B.rows <= 0 || B.cols <= 0) {
    return false;
  }
  if (A.col_ptrs.size() != static_cast<size_t>(A.cols) + 1 || B.col_ptrs.size() != static_cast<size_t>(B.cols) + 1) {
    return false;
  }
  if (A.col_ptrs.empty() || A.col_ptrs[0] != 0 || B.col_ptrs.empty() || B.col_ptrs[0] != 0) {
    return false;
  }
  if (std::cmp_not_equal(A.nnz, A.values.size()) || std::cmp_not_equal(B.nnz, B.values.size())) {
    return false;
  }
  return true;
}

bool BarkalovaMMultMatrixCcsSTL::PreProcessingImpl() {
  return true;
}

namespace {
constexpr double kEpsilon = 1e-10;

void TransponirMatr(const CCSMatrix &a, CCSMatrix &at) {
  at.rows = a.cols;
  at.cols = a.rows;
  at.nnz = a.nnz;

  if (a.nnz == 0) {
    at.values.clear();
    at.row_indices.clear();
    at.col_ptrs.assign(at.cols + 1, 0);
    return;
  }

  std::vector<int> row_count(at.cols, 0);
  for (int i = 0; i < a.nnz; i++) {
    row_count[a.row_indices[i]]++;
  }

  at.col_ptrs.resize(at.cols + 1);
  at.col_ptrs[0] = 0;
  for (int i = 0; i < at.cols; i++) {
    at.col_ptrs[i + 1] = at.col_ptrs[i] + row_count[i];
  }

  at.values.resize(a.nnz);
  at.row_indices.resize(a.nnz);

  std::vector<int> current_pos(at.cols, 0);
  for (int col = 0; col < a.cols; col++) {
    for (int i = a.col_ptrs[col]; i < a.col_ptrs[col + 1]; i++) {
      int row = a.row_indices[i];
      Complex val = a.values[i];

      int pos = at.col_ptrs[row] + current_pos[row];
      at.values[pos] = val;
      at.row_indices[pos] = col;
      current_pos[row]++;
    }
  }
}

bool IsNonZero(const Complex &val) {
  return std::abs(val.real()) > kEpsilon || std::abs(val.imag()) > kEpsilon;
}

Complex ComputeScalarProduct(const CCSMatrix &at, const CCSMatrix &b, int row_a, int col_b) {
  Complex sum = Complex(0.0, 0.0);

  int ks = at.col_ptrs[row_a];
  int ls = b.col_ptrs[col_b];
  int kf = at.col_ptrs[row_a + 1];
  int lf = b.col_ptrs[col_b + 1];

  while ((ks < kf) && (ls < lf)) {
    if (at.row_indices[ks] < b.row_indices[ls]) {
      ks++;
    } else if (at.row_indices[ks] > b.row_indices[ls]) {
      ls++;
    } else {
      sum += at.values[ks] * b.values[ls];
      ks++;
      ls++;
    }
  }

  return sum;
}

struct ColumnResult {
  std::vector<int> rows;
  std::vector<Complex> values;
};

void ProcessColumnRange(int start_col, int end_col, const CCSMatrix &at, const CCSMatrix &b,
                        std::vector<std::vector<int>> &col_rows, std::vector<std::vector<Complex>> &col_vals) {
  for (int col = start_col; col < end_col; ++col) {
    std::vector<int> rows;
    std::vector<Complex> vals;
    rows.reserve(100);
    vals.reserve(100);

    for (int i = 0; i < at.cols; i++) {
      Complex sum = ComputeScalarProduct(at, b, i, col);
      if (IsNonZero(sum)) {
        rows.push_back(i);
        vals.push_back(sum);
      }
    }

    col_rows[col] = std::move(rows);
    col_vals[col] = std::move(vals);
  }
}

}  // namespace

bool BarkalovaMMultMatrixCcsSTL::RunImpl() {
  const auto &a = GetInput().first;
  const auto &b = GetInput().second;

  try {
    CCSMatrix at;
    TransponirMatr(a, at);
    CCSMatrix c;
    c.rows = a.rows;
    c.cols = b.cols;

    const int total_cols = c.cols;
    unsigned int hardware_threads = std::thread::hardware_concurrency();
    const unsigned int num_threads = std::max(1U, hardware_threads);

    std::vector<std::vector<int>> col_rows(total_cols);
    std::vector<std::vector<Complex>> col_vals(total_cols);
    std::vector<std::future<void>> futures;
    futures.reserve(num_threads);

    int cols_per_thread = total_cols / static_cast<int>(num_threads);
    int remainder = total_cols % static_cast<int>(num_threads);
    int current_start = 0;

    for (unsigned int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
      int cols_for_thread = cols_per_thread + (std::cmp_less(thread_idx, remainder) ? 1 : 0);
      int start = current_start;
      int end = current_start + cols_for_thread;
      current_start = end;
      if (start >= total_cols) {
        break;
      }

      futures.push_back(std::async(std::launch::async, ProcessColumnRange, start, end, std::cref(at), std::cref(b),
                                   std::ref(col_rows), std::ref(col_vals)));
    }
    for (auto &future : futures) {
      future.get();
    }

    std::vector<int> col_ptrs = {0};
    std::vector<int> row_indices;
    std::vector<Complex> values;

    for (int j = 0; j < total_cols; ++j) {
      row_indices.insert(row_indices.end(), col_rows[j].begin(), col_rows[j].end());
      values.insert(values.end(), col_vals[j].begin(), col_vals[j].end());
      col_ptrs.push_back(static_cast<int>(values.size()));
    }

    c.values = std::move(values);
    c.row_indices = std::move(row_indices);
    c.col_ptrs = std::move(col_ptrs);
    c.nnz = static_cast<int>(c.values.size());

    GetOutput() = c;
    return true;

  } catch (const std::exception &) {
    return false;
  }
}

bool BarkalovaMMultMatrixCcsSTL::PostProcessingImpl() {
  const auto &c = GetOutput();
  return c.rows > 0 && c.cols > 0 && c.col_ptrs.size() == static_cast<size_t>(c.cols) + 1;
}

}  // namespace barkalova_m_mult_matrix_ccs
