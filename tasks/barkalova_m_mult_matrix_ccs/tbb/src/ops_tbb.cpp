#include "barkalova_m_mult_matrix_ccs/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include <cmath>
#include <complex>
#include <cstddef>
#include <exception>
#include <utility>
#include <vector>

#include "barkalova_m_mult_matrix_ccs/common/include/common.hpp"

namespace barkalova_m_mult_matrix_ccs {

BarkalovaMMultMatrixCcsTBB::BarkalovaMMultMatrixCcsTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = CCSMatrix{};
}

bool BarkalovaMMultMatrixCcsTBB::ValidationImpl() {
  const auto &[A, B] = GetInput();
  if (A.cols != B.rows) {
    return false;
  }
  if (A.rows <= 0 || A.cols <= 0 || B.rows <= 0 || B.cols <= 0) {
    return false;
  }
  if (A.col_ptrs.size() != static_cast<size_t>(A.cols) + 1) {
    return false;
  }
  if (B.col_ptrs.size() != static_cast<size_t>(B.cols) + 1) {
    return false;
  }
  if (A.col_ptrs[0] != 0 || B.col_ptrs[0] != 0) {
    return false;
  }
  if (static_cast<size_t>(A.nnz) != A.values.size()) {
    return false;
  }
  if (static_cast<size_t>(B.nnz) != B.values.size()) {
    return false;
  }
  return true;
}

bool BarkalovaMMultMatrixCcsTBB::PreProcessingImpl() {
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

void ProcessColumn(int j, const CCSMatrix &at, const CCSMatrix &b, std::vector<int> &out_rows,
                   std::vector<Complex> &out_vals) {
  out_rows.reserve(100);
  out_vals.reserve(100);

  for (int i = 0; i < at.cols; i++) {
    Complex sum = Complex(0.0, 0.0);

    int ks = at.col_ptrs[i];
    int ls = b.col_ptrs[j];
    int kf = at.col_ptrs[i + 1];
    int lf = b.col_ptrs[j + 1];

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

    if (IsNonZero(sum)) {
      out_rows.push_back(i);
      out_vals.push_back(sum);
    }
  }
}

}  // namespace

bool BarkalovaMMultMatrixCcsTBB::RunImpl() {
  const auto &a = GetInput().first;
  const auto &b = GetInput().second;

  try {
    CCSMatrix at;
    TransponirMatr(a, at);

    CCSMatrix c;
    c.rows = a.rows;
    c.cols = b.cols;

    std::vector<std::vector<int>> col_rows(c.cols);
    std::vector<std::vector<Complex>> col_vals(c.cols);

    tbb::parallel_for(tbb::blocked_range<int>(0, c.cols), [&](const tbb::blocked_range<int> &range) {
      for (int j = range.begin(); j < range.end(); ++j) {
        ProcessColumn(j, at, b, col_rows[j], col_vals[j]);
      }
    });

    std::vector<int> col_ptrs = {0};
    std::vector<int> row_indices;
    std::vector<Complex> values;

    for (int j = 0; j < c.cols; j++) {
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

bool BarkalovaMMultMatrixCcsTBB::PostProcessingImpl() {
  const auto &c = GetOutput();
  return c.rows > 0 && c.cols > 0 && c.col_ptrs.size() == static_cast<size_t>(c.cols) + 1;
}

}  // namespace barkalova_m_mult_matrix_ccs
