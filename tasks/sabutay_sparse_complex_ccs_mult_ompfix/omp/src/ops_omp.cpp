#include "sabutay_sparse_complex_ccs_mult_ompfix/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <utility>
#include <vector>

#include "sabutay_sparse_complex_ccs_mult_ompfix/common/include/common.hpp"

namespace sabutay_sparse_complex_ccs_mult_ompfix {
namespace {

constexpr double kDropMagnitude = 1e-12;

auto IsValidStructure(const CCS &matrix) -> bool {
  if (matrix.row_count < 0 || matrix.col_count < 0) {
    return false;
  }
  if (matrix.col_start.size() != (static_cast<std::size_t>(matrix.col_count) + 1U)) {
    return false;
  }
  if (matrix.row_index.size() != matrix.nz.size()) {
    return false;
  }
  if (matrix.col_start.empty() || matrix.col_start.front() != 0) {
    return false;
  }
  if (!std::cmp_equal(matrix.col_start.back(), matrix.nz.size())) {
    return false;
  }
  for (int j = 0; j < matrix.col_count; ++j) {
    const auto col_idx = static_cast<std::size_t>(j);
    if (matrix.col_start[col_idx] > matrix.col_start[col_idx + 1U]) {
      return false;
    }
  }
  return std::ranges::all_of(matrix.row_index, [&matrix](int row) { return row >= 0 && row < matrix.row_count; });
}

void BuildColumnFromRight(const CCS &left, const CCS &right, int column_index,
                          std::vector<std::pair<int, std::complex<double>>> &buffer) {
  const int right_begin = right.col_start[static_cast<std::size_t>(column_index)];
  const int right_end = right.col_start[static_cast<std::size_t>(column_index) + 1U];
  buffer.clear();

  for (int right_pos = right_begin; right_pos < right_end; ++right_pos) {
    const int inner = right.row_index[static_cast<std::size_t>(right_pos)];
    const std::complex<double> scalar = right.nz[static_cast<std::size_t>(right_pos)];
    const int left_begin = left.col_start[static_cast<std::size_t>(inner)];
    const int left_end = left.col_start[static_cast<std::size_t>(inner) + 1U];
    for (int left_pos = left_begin; left_pos < left_end; ++left_pos) {
      const int row = left.row_index[static_cast<std::size_t>(left_pos)];
      buffer.emplace_back(row, left.nz[static_cast<std::size_t>(left_pos)] * scalar);
    }
  }
}

void CoalesceBufferToColumn(const std::vector<std::pair<int, std::complex<double>>> &buffer,
                            std::vector<std::pair<int, std::complex<double>>> &column) {
  if (buffer.empty()) {
    column.clear();
    return;
  }

  column.clear();
  column.reserve(buffer.size());
  int active_row = buffer[0].first;
  std::complex<double> running = buffer[0].second;
  for (std::size_t idx = 1; idx < buffer.size(); ++idx) {
    const auto &entry = buffer[idx];
    if (entry.first == active_row) {
      running += entry.second;
    } else {
      if (std::abs(running) > kDropMagnitude) {
        column.emplace_back(active_row, running);
      }
      active_row = entry.first;
      running = entry.second;
    }
  }
  if (std::abs(running) > kDropMagnitude) {
    column.emplace_back(active_row, running);
  }
}

void SortBufferByRow(std::vector<std::pair<int, std::complex<double>>> &buffer) {
  std::ranges::sort(buffer, {}, &std::pair<int, std::complex<double>>::first);
}

}  // namespace

SabutaySparseComplexCcsMultOmpFix::SabutaySparseComplexCcsMultOmpFix(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = CCS();
}

bool SabutaySparseComplexCcsMultOmpFix::ValidationImpl() {
  const CCS &left = std::get<0>(GetInput());
  const CCS &right = std::get<1>(GetInput());
  return left.col_count == right.row_count && IsValidStructure(left) && IsValidStructure(right);
}

bool SabutaySparseComplexCcsMultOmpFix::PreProcessingImpl() {
  return true;
}

void SabutaySparseComplexCcsMultOmpFix::BuildProductMatrix(const CCS &left, const CCS &right, CCS &out) {
  out.row_count = left.row_count;
  out.col_count = right.col_count;
  out.col_start.assign(static_cast<std::size_t>(out.col_count) + 1U, 0);
  out.row_index.clear();
  out.nz.clear();
  if (out.col_count == 0) {
    return;
  }

  std::vector<std::vector<std::pair<int, std::complex<double>>>> columns(static_cast<std::size_t>(out.col_count));
#pragma omp parallel default(none) shared(left, right, columns)
  {
    std::vector<std::pair<int, std::complex<double>>> buffer;
    buffer.reserve(128U);

#pragma omp for schedule(static)
    for (int j = 0; j < right.col_count; ++j) {
      auto &column = columns[static_cast<std::size_t>(j)];
      BuildColumnFromRight(left, right, j, buffer);
      if (!buffer.empty()) {
        SortBufferByRow(buffer);
      }
      CoalesceBufferToColumn(buffer, column);
    }
  }

  for (int j = 0; j < out.col_count; ++j) {
    const int col_size = static_cast<int>(columns[static_cast<std::size_t>(j)].size());
    out.col_start[static_cast<std::size_t>(j) + 1U] = out.col_start[static_cast<std::size_t>(j)] + col_size;
  }

  const auto nnz = static_cast<std::size_t>(out.col_start.back());
  out.row_index.resize(nnz);
  out.nz.resize(nnz);

  for (int j = 0; j < out.col_count; ++j) {
    const int start = out.col_start[static_cast<std::size_t>(j)];
    const auto &column = columns[static_cast<std::size_t>(j)];
    for (std::size_t k = 0; k < column.size(); ++k) {
      const int dst = start + static_cast<int>(k);
      out.row_index[static_cast<std::size_t>(dst)] = column[k].first;
      out.nz[static_cast<std::size_t>(dst)] = column[k].second;
    }
  }
}

bool SabutaySparseComplexCcsMultOmpFix::RunImpl() {
  const CCS &left = std::get<0>(GetInput());
  const CCS &right = std::get<1>(GetInput());
  BuildProductMatrix(left, right, GetOutput());
  return true;
}

bool SabutaySparseComplexCcsMultOmpFix::PostProcessingImpl() {
  return true;
}

}  // namespace sabutay_sparse_complex_ccs_mult_ompfix
