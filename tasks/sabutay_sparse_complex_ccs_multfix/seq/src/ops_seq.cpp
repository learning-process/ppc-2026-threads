#include "sabutay_sparse_complex_ccs_multfix/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <utility>
#include <vector>

#include "sabutay_sparse_complex_ccs_multfix/common/include/common.hpp"

namespace sabutay_sparse_complex_ccs_multfix {
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

void SpmmAbc(const CCS &a, const CCS &b, CCS &c) {
  c.row_count = a.row_count;
  c.col_count = b.col_count;
  c.col_start.assign(static_cast<std::size_t>(c.col_count) + 1U, 0);
  c.row_index.clear();
  c.nz.clear();
  if (c.col_count == 0) {
    return;
  }

  std::vector<std::pair<int, std::complex<double>>> buffer;
  buffer.reserve(128U);

  for (int j = 0; j < b.col_count; ++j) {
    const int b_begin = b.col_start[static_cast<std::size_t>(j)];
    const int b_end = b.col_start[static_cast<std::size_t>(j) + 1U];
    buffer.clear();
    for (int b_pos = b_begin; b_pos < b_end; ++b_pos) {
      const int k = b.row_index[static_cast<std::size_t>(b_pos)];
      const std::complex<double> s = b.nz[static_cast<std::size_t>(b_pos)];
      const int a_lo = a.col_start[static_cast<std::size_t>(k)];
      const int a_hi = a.col_start[static_cast<std::size_t>(k) + 1U];
      for (int s_idx = a_lo; s_idx < a_hi; ++s_idx) {
        const int i = a.row_index[static_cast<std::size_t>(s_idx)];
        buffer.emplace_back(i, a.nz[static_cast<std::size_t>(s_idx)] * s);
      }
    }
    if (buffer.empty()) {
      c.col_start[static_cast<std::size_t>(j) + 1U] = static_cast<int>(c.nz.size());
      continue;
    }
    std::ranges::sort(buffer, {}, &std::pair<int, std::complex<double>>::first);
    CoalesceSortedPairs(buffer, c);
    c.col_start[static_cast<std::size_t>(j) + 1U] = static_cast<int>(c.nz.size());
  }
}

}  // namespace

SabutaySparseComplexCcsMultFixSEQ::SabutaySparseComplexCcsMultFixSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = CCS();
}

void SabutaySparseComplexCcsMultFixSEQ::BuildProductMatrix(const CCS &left, const CCS &right, CCS &out) {
  SpmmAbc(left, right, out);
}

bool SabutaySparseComplexCcsMultFixSEQ::ValidationImpl() {
  const CCS &left = std::get<0>(GetInput());
  const CCS &right = std::get<1>(GetInput());
  return left.col_count == right.row_count;
}

bool SabutaySparseComplexCcsMultFixSEQ::PreProcessingImpl() {
  return true;
}

bool SabutaySparseComplexCcsMultFixSEQ::RunImpl() {
  const CCS &left = std::get<0>(GetInput());
  const CCS &right = std::get<1>(GetInput());
  BuildProductMatrix(left, right, GetOutput());
  return true;
}

bool SabutaySparseComplexCcsMultFixSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace sabutay_sparse_complex_ccs_multfix
