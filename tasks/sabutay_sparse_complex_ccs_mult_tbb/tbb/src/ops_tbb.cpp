#include "sabutay_sparse_complex_ccs_mult_tbb/tbb/include/ops_tbb.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <utility>
#include <vector>

#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/parallel_for.h"
#include "sabutay_sparse_complex_ccs_mult_tbb/common/include/common.hpp"

namespace sabutay_sparse_complex_ccs_mult_tbb {

namespace {

constexpr double kDropMagnitude = 1e-14;

void CoalesceSortedPairs(const std::vector<std::pair<int, std::complex<double>>> &row_sorted, CCS &out) {
  if (row_sorted.empty()) {
    return;
  }
  int active_row = row_sorted[0].first;
  std::complex<double> running = row_sorted[0].second;
  for (std::size_t idx = 1; idx < row_sorted.size(); ++idx) {
    const int r = row_sorted[idx].first;
    if (r == active_row) {
      running += row_sorted[idx].second;
    } else {
      if (std::abs(running) > kDropMagnitude) {
        out.row_index.push_back(active_row);
        out.nz.push_back(running);
      }
      active_row = r;
      running = row_sorted[idx].second;
    }
  }
  if (std::abs(running) > kDropMagnitude) {
    out.row_index.push_back(active_row);
    out.nz.push_back(running);
  }
}

void BuildColumn(const CCS &left, const CCS &right, int jcol,
                 std::vector<std::pair<int, std::complex<double>>> &buffer) {
  const int b_begin = right.col_start[static_cast<std::size_t>(jcol)];
  const int b_end = right.col_start[static_cast<std::size_t>(jcol) + 1U];
  buffer.clear();
  for (int b_pos = b_begin; b_pos < b_end; ++b_pos) {
    const int k = right.row_index[static_cast<std::size_t>(b_pos)];
    const std::complex<double> s = right.nz[static_cast<std::size_t>(b_pos)];
    const int a_lo = left.col_start[static_cast<std::size_t>(k)];
    const int a_hi = left.col_start[static_cast<std::size_t>(k) + 1U];
    for (int s_idx = a_lo; s_idx < a_hi; ++s_idx) {
      const int i = left.row_index[static_cast<std::size_t>(s_idx)];
      buffer.emplace_back(i, left.nz[static_cast<std::size_t>(s_idx)] * s);
    }
  }
}

bool IsValidCCS(const CCS &matrix) {
  if (matrix.row_count < 0 || matrix.col_count < 0) {
    return false;
  }
  if (matrix.col_start.size() != static_cast<std::size_t>(matrix.col_count) + 1U) {
    return false;
  }
  if (matrix.row_index.size() != matrix.nz.size()) {
    return false;
  }
  if (matrix.col_start.empty() || matrix.col_start[0] != 0) {
    return false;
  }
  if (static_cast<std::size_t>(matrix.col_start.back()) != matrix.row_index.size()) {
    return false;
  }
  for (int jcol = 0; jcol < matrix.col_count; ++jcol) {
    if (matrix.col_start[static_cast<std::size_t>(jcol)] > matrix.col_start[static_cast<std::size_t>(jcol) + 1U]) {
      return false;
    }
  }
  return std::ranges::all_of(matrix.row_index, [&](int row) { return row >= 0 && row < matrix.row_count; });
}

}  // namespace

SabutaySparseComplexCcsMultFixTBB::SabutaySparseComplexCcsMultFixTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = CCS();
}

void SabutaySparseComplexCcsMultFixTBB::BuildProductMatrix(const CCS &left, const CCS &right, CCS &out) {
  out.row_count = left.row_count;
  out.col_count = right.col_count;
  out.col_start.assign(static_cast<std::size_t>(out.col_count) + 1U, 0);
  out.row_index.clear();
  out.nz.clear();
  if (out.col_count == 0) {
    return;
  }

  std::vector<std::vector<int>> local_row_index(static_cast<std::size_t>(right.col_count));
  std::vector<std::vector<std::complex<double>>> local_nz(static_cast<std::size_t>(right.col_count));
  std::vector<int> local_sizes(static_cast<std::size_t>(right.col_count), 0);

  oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<int>(0, right.col_count),
                            [&](const oneapi::tbb::blocked_range<int> &range) {
    std::vector<std::pair<int, std::complex<double>>> buffer;
    buffer.reserve(128U);
    for (int jcol = range.begin(); jcol < range.end(); ++jcol) {
      BuildColumn(left, right, jcol, buffer);
      if (!buffer.empty()) {
        std::ranges::sort(buffer, {}, &std::pair<int, std::complex<double>>::first);
        CCS tmp;
        CoalesceSortedPairs(buffer, tmp);
        local_row_index[static_cast<std::size_t>(jcol)] = std::move(tmp.row_index);
        local_nz[static_cast<std::size_t>(jcol)] = std::move(tmp.nz);
      } else {
        local_row_index[static_cast<std::size_t>(jcol)].clear();
        local_nz[static_cast<std::size_t>(jcol)].clear();
      }
      local_sizes[static_cast<std::size_t>(jcol)] = static_cast<int>(local_nz[static_cast<std::size_t>(jcol)].size());
    }
  });

  for (int jcol = 0; jcol < right.col_count; ++jcol) {
    out.col_start[static_cast<std::size_t>(jcol) + 1U] =
        out.col_start[static_cast<std::size_t>(jcol)] + local_sizes[static_cast<std::size_t>(jcol)];
  }

  out.row_index.reserve(static_cast<std::size_t>(out.col_start[static_cast<std::size_t>(out.col_count)]));
  out.nz.reserve(static_cast<std::size_t>(out.col_start[static_cast<std::size_t>(out.col_count)]));

  for (int jcol = 0; jcol < right.col_count; ++jcol) {
    const auto idx = static_cast<std::size_t>(jcol);
    out.row_index.insert(out.row_index.end(), local_row_index[idx].begin(), local_row_index[idx].end());
    out.nz.insert(out.nz.end(), local_nz[idx].begin(), local_nz[idx].end());
  }
}

bool SabutaySparseComplexCcsMultFixTBB::ValidationImpl() {
  const CCS &left = std::get<0>(GetInput());
  const CCS &right = std::get<1>(GetInput());
  return IsValidCCS(left) && IsValidCCS(right) && left.col_count == right.row_count;
}

bool SabutaySparseComplexCcsMultFixTBB::PreProcessingImpl() {
  return true;
}

bool SabutaySparseComplexCcsMultFixTBB::RunImpl() {
  const CCS &left = std::get<0>(GetInput());
  const CCS &right = std::get<1>(GetInput());
  BuildProductMatrix(left, right, GetOutput());
  return true;
}

bool SabutaySparseComplexCcsMultFixTBB::PostProcessingImpl() {
  return true;
}

}  // namespace sabutay_sparse_complex_ccs_mult_tbb
