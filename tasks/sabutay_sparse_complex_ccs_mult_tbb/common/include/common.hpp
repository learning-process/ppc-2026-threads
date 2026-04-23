#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "task/include/task.hpp"

namespace sabutay_sparse_complex_ccs_mult_tbb {

using Complex = std::complex<double>;
using MatrixEntry = std::pair<int, Complex>;

struct SparseMatrixCCS {
  int rows = 0;
  int cols = 0;
  std::vector<int> col_ptr = {0};
  std::vector<int> row_ind;
  std::vector<Complex> values;
};

inline constexpr double kComplexEps = 1e-12;

[[nodiscard]] inline bool IsNearZero(const Complex &value, double eps = kComplexEps) {
  return std::norm(value) <= (eps * eps);
}

[[nodiscard]] inline SparseMatrixCCS MakeZeroMatrix(int rows, int cols) {
  if (rows < 0 || cols < 0) {
    return {};
  }

  SparseMatrixCCS matrix;
  matrix.rows = rows;
  matrix.cols = cols;
  matrix.col_ptr.assign(static_cast<std::size_t>(cols) + 1, 0);
  return matrix;
}

[[nodiscard]] inline bool HasValidCcsShape(const SparseMatrixCCS &matrix) {
  return matrix.rows >= 0 && matrix.cols >= 0 &&
         matrix.col_ptr.size() == (static_cast<std::size_t>(matrix.cols) + 1U) &&
         matrix.row_ind.size() == matrix.values.size() && !matrix.col_ptr.empty() && matrix.col_ptr.front() == 0;
}

[[nodiscard]] inline bool HasValidColumnPointers(const SparseMatrixCCS &matrix, int nnz) {
  if (matrix.col_ptr.back() != nnz) {
    return false;
  }

  for (std::size_t col = 0; (col + 1U) < matrix.col_ptr.size(); ++col) {
    const int begin = matrix.col_ptr[col];
    const int end = matrix.col_ptr[col + 1U];
    if (begin > end || begin < 0 || end > nnz) {
      return false;
    }
  }

  return true;
}

[[nodiscard]] inline bool HasValidRowIndices(const SparseMatrixCCS &matrix) {
  return std::ranges::all_of(matrix.row_ind, [&matrix](const int row) { return row >= 0 && row < matrix.rows; });
}

[[nodiscard]] inline bool IsValidCcs(const SparseMatrixCCS &matrix) {
  if (!HasValidCcsShape(matrix)) {
    return false;
  }

  const auto nnz = static_cast<int>(matrix.row_ind.size());
  return HasValidColumnPointers(matrix, nnz) && HasValidRowIndices(matrix);
}

inline void SortEntriesByRow(std::vector<MatrixEntry> &entries) {
  std::ranges::sort(entries, {}, &MatrixEntry::first);
}

inline void FlushPendingEntry(std::optional<MatrixEntry> &pending, std::vector<int> &row_ind,
                              std::vector<Complex> &values, double eps) {
  if (pending.has_value() && !IsNearZero(pending->second, eps)) {
    row_ind.push_back(pending->first);
    values.push_back(pending->second);
  }
}

inline void AppendCompressedEntries(const std::vector<MatrixEntry> &entries, std::vector<int> &row_ind,
                                    std::vector<Complex> &values, double eps) {
  std::optional<MatrixEntry> pending;

  for (const auto &[row, value] : entries) {
    if (pending.has_value() && pending->first == row) {
      pending->second += value;
      continue;
    }

    FlushPendingEntry(pending, row_ind, values, eps);
    pending = MatrixEntry{row, value};
  }

  FlushPendingEntry(pending, row_ind, values, eps);
}

[[nodiscard]] inline std::vector<MatrixEntry> GetColumnEntries(const SparseMatrixCCS &matrix, int col) {
  const int begin = matrix.col_ptr[static_cast<std::size_t>(col)];
  const int end = matrix.col_ptr[static_cast<std::size_t>(col) + 1U];
  std::vector<MatrixEntry> entries;
  entries.reserve(static_cast<std::size_t>(end - begin));

  for (int idx = begin; idx < end; ++idx) {
    entries.emplace_back(matrix.row_ind[static_cast<std::size_t>(idx)], matrix.values[static_cast<std::size_t>(idx)]);
  }

  SortEntriesByRow(entries);
  return entries;
}

[[nodiscard]] inline SparseMatrixCCS NormalizeCcs(const SparseMatrixCCS &matrix, double eps = kComplexEps) {
  if (!IsValidCcs(matrix)) {
    return {};
  }

  SparseMatrixCCS normalized = MakeZeroMatrix(matrix.rows, matrix.cols);
  normalized.values.reserve(matrix.values.size());
  normalized.row_ind.reserve(matrix.row_ind.size());

  for (int col = 0; col < matrix.cols; ++col) {
    const auto entries = GetColumnEntries(matrix, col);
    AppendCompressedEntries(entries, normalized.row_ind, normalized.values, eps);
    normalized.col_ptr[static_cast<std::size_t>(col) + 1U] = static_cast<int>(normalized.row_ind.size());
  }

  return normalized;
}

[[nodiscard]] inline std::vector<std::vector<Complex>> ToDense(const SparseMatrixCCS &matrix) {
  if (!IsValidCcs(matrix)) {
    return {};
  }

  std::vector<std::vector<Complex>> dense(static_cast<std::size_t>(matrix.rows),
                                          std::vector<Complex>(static_cast<std::size_t>(matrix.cols), Complex{}));

  for (int col = 0; col < matrix.cols; ++col) {
    const int begin = matrix.col_ptr[static_cast<std::size_t>(col)];
    const int end = matrix.col_ptr[static_cast<std::size_t>(col) + 1U];

    for (int idx = begin; idx < end; ++idx) {
      dense[static_cast<std::size_t>(matrix.row_ind[static_cast<std::size_t>(idx)])][static_cast<std::size_t>(col)] =
          matrix.values[static_cast<std::size_t>(idx)];
    }
  }

  return dense;
}

[[nodiscard]] inline SparseMatrixCCS DenseToCcs(const std::vector<std::vector<Complex>> &dense,
                                                double eps = kComplexEps) {
  const int rows = static_cast<int>(dense.size());
  const int cols = rows == 0 ? 0 : static_cast<int>(dense.front().size());

  for (const auto &row : dense) {
    if (std::cmp_not_equal(row.size(), static_cast<std::size_t>(cols))) {
      return {};
    }
  }

  SparseMatrixCCS matrix = MakeZeroMatrix(rows, cols);

  for (int col = 0; col < cols; ++col) {
    for (int row = 0; row < rows; ++row) {
      const Complex &value = dense[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)];
      if (!IsNearZero(value, eps)) {
        matrix.row_ind.push_back(row);
        matrix.values.push_back(value);
      }
    }
    matrix.col_ptr[static_cast<std::size_t>(col) + 1U] = static_cast<int>(matrix.row_ind.size());
  }

  return matrix;
}

inline void AccumulateProductRows(const SparseMatrixCCS &lhs, const SparseMatrixCCS &rhs, int col,
                                  std::vector<Complex> &accumulator, std::vector<int> &marker,
                                  std::vector<int> &touched_rows) {
  for (int rhs_idx = rhs.col_ptr[static_cast<std::size_t>(col)];
       rhs_idx < rhs.col_ptr[static_cast<std::size_t>(col) + 1U]; ++rhs_idx) {
    const int lhs_col = rhs.row_ind[static_cast<std::size_t>(rhs_idx)];
    const Complex rhs_value = rhs.values[static_cast<std::size_t>(rhs_idx)];

    for (int lhs_idx = lhs.col_ptr[static_cast<std::size_t>(lhs_col)];
         lhs_idx < lhs.col_ptr[static_cast<std::size_t>(lhs_col) + 1U]; ++lhs_idx) {
      const int row = lhs.row_ind[static_cast<std::size_t>(lhs_idx)];

      if (marker[static_cast<std::size_t>(row)] != col) {
        marker[static_cast<std::size_t>(row)] = col;
        accumulator[static_cast<std::size_t>(row)] = Complex{};
        touched_rows.push_back(row);
      }

      accumulator[static_cast<std::size_t>(row)] += lhs.values[static_cast<std::size_t>(lhs_idx)] * rhs_value;
    }
  }
}

inline void AppendAccumulatedRows(const std::vector<Complex> &accumulator, std::vector<int> &touched_rows, double eps,
                                  SparseMatrixCCS &result) {
  std::ranges::sort(touched_rows);

  for (const int row : touched_rows) {
    const Complex value = accumulator[static_cast<std::size_t>(row)];
    if (!IsNearZero(value, eps)) {
      result.row_ind.push_back(row);
      result.values.push_back(value);
    }
  }
}

[[nodiscard]] inline SparseMatrixCCS MultiplyCcsReference(const SparseMatrixCCS &lhs, const SparseMatrixCCS &rhs,
                                                          double eps = kComplexEps) {
  if (!IsValidCcs(lhs) || !IsValidCcs(rhs) || lhs.cols != rhs.rows) {
    return {};
  }

  const SparseMatrixCCS left = NormalizeCcs(lhs, eps);
  const SparseMatrixCCS right = NormalizeCcs(rhs, eps);
  SparseMatrixCCS result = MakeZeroMatrix(left.rows, right.cols);

  std::vector<Complex> accumulator(static_cast<std::size_t>(left.rows), Complex{});
  std::vector<int> marker(static_cast<std::size_t>(left.rows), -1);
  std::vector<int> touched_rows;

  for (int col = 0; col < right.cols; ++col) {
    touched_rows.clear();
    AccumulateProductRows(left, right, col, accumulator, marker, touched_rows);
    AppendAccumulatedRows(accumulator, touched_rows, eps, result);
    result.col_ptr[static_cast<std::size_t>(col) + 1U] = static_cast<int>(result.row_ind.size());
  }

  return result;
}

[[nodiscard]] inline bool AreMatricesEqual(const SparseMatrixCCS &lhs, const SparseMatrixCCS &rhs,
                                           double eps = kComplexEps) {
  if (!IsValidCcs(lhs) || !IsValidCcs(rhs)) {
    return false;
  }

  const SparseMatrixCCS left = NormalizeCcs(lhs, eps);
  const SparseMatrixCCS right = NormalizeCcs(rhs, eps);

  if (left.rows != right.rows || left.cols != right.cols) {
    return false;
  }

  if (left.col_ptr != right.col_ptr || left.row_ind != right.row_ind) {
    return false;
  }

  if (left.values.size() != right.values.size()) {
    return false;
  }

  for (std::size_t idx = 0; idx < left.values.size(); ++idx) {
    if (!IsNearZero(left.values[idx] - right.values[idx], eps)) {
      return false;
    }
  }

  return true;
}

[[nodiscard]] inline SparseMatrixCCS BuildDeterministicMatrix(int rows, int cols, int non_zero_per_col,
                                                              int offset = 0) {
  if (rows < 0 || cols < 0 || non_zero_per_col < 0) {
    return {};
  }

  SparseMatrixCCS matrix = MakeZeroMatrix(rows, cols);
  if (rows == 0 || cols == 0 || non_zero_per_col == 0) {
    return matrix;
  }

  for (int col = 0; col < cols; ++col) {
    const int count = std::min(rows, non_zero_per_col);
    std::vector<MatrixEntry> entries;
    entries.reserve(static_cast<std::size_t>(count));

    for (int idx = 0; idx < count; ++idx) {
      const int row = (offset + (col * 7) + (idx * 11)) % rows;
      const auto real = static_cast<double>((((col + 1) * (idx + 2 + offset)) % 17) - 8);
      const auto imag = static_cast<double>((((col + 3 + offset) * (idx + 1)) % 13) - 6);
      entries.emplace_back(row, Complex(real, imag));
    }

    SortEntriesByRow(entries);
    AppendCompressedEntries(entries, matrix.row_ind, matrix.values, kComplexEps);
    matrix.col_ptr[static_cast<std::size_t>(col) + 1U] = static_cast<int>(matrix.row_ind.size());
  }

  return matrix;
}

using InType = std::tuple<SparseMatrixCCS, SparseMatrixCCS>;
using OutType = SparseMatrixCCS;
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace sabutay_sparse_complex_ccs_mult_tbb
