#include "sabutay_sparse_complex_ccs_mult_all/all/include/ops_all.hpp"

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <complex>
#include <cstddef>
#include <thread>
#include <utility>
#include <vector>

#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/parallel_for.h"
#include "sabutay_sparse_complex_ccs_mult_all/common/include/common.hpp"
#include "task/include/task.hpp"
#include "util/include/util.hpp"

namespace sabutay_sparse_complex_ccs_mult_all {
namespace {

constexpr double kDropMagnitude = 1e-14;

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

void SortBufferByRow(std::vector<std::pair<int, std::complex<double>>> &buffer) {
  std::ranges::sort(buffer, {}, &std::pair<int, std::complex<double>>::first);
}

void CoalesceSortedPairs(const std::vector<std::pair<int, std::complex<double>>> &row_sorted, CCS &out,
                         double drop_magnitude) {
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
      if (std::abs(running) > drop_magnitude) {
        out.row_index.push_back(active_row);
        out.nz.push_back(running);
      }
      active_row = r;
      running = row_sorted[idx].second;
    }
  }
  if (std::abs(running) > drop_magnitude) {
    out.row_index.push_back(active_row);
    out.nz.push_back(running);
  }
}

void CoalesceBufferToColumn(const std::vector<std::pair<int, std::complex<double>>> &buffer,
                            std::vector<std::pair<int, std::complex<double>>> &column, double drop_magnitude) {
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
      if (std::abs(running) > drop_magnitude) {
        column.emplace_back(active_row, running);
      }
      active_row = entry.first;
      running = entry.second;
    }
  }
  if (std::abs(running) > drop_magnitude) {
    column.emplace_back(active_row, running);
  }
}

void SpmmAbc(const CCS &left, const CCS &right, CCS &out, double drop_magnitude) {
  out.row_count = left.row_count;
  out.col_count = right.col_count;
  out.col_start.assign(static_cast<std::size_t>(out.col_count) + 1U, 0);
  out.row_index.clear();
  out.nz.clear();
  if (out.col_count == 0) {
    return;
  }

  std::vector<std::pair<int, std::complex<double>>> buffer;
  buffer.reserve(128U);

  for (int j = 0; j < right.col_count; ++j) {
    BuildColumnFromRight(left, right, j, buffer);
    if (buffer.empty()) {
      out.col_start[static_cast<std::size_t>(j) + 1U] = static_cast<int>(out.nz.size());
      continue;
    }
    SortBufferByRow(buffer);
    CoalesceSortedPairs(buffer, out, drop_magnitude);
    out.col_start[static_cast<std::size_t>(j) + 1U] = static_cast<int>(out.nz.size());
  }
}

void BuildProductMatrixOmp(const CCS &left, const CCS &right, CCS &out) {
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
      CoalesceBufferToColumn(buffer, column, kDropMagnitude);
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

void BuildProductMatrixTbb(const CCS &left, const CCS &right, CCS &out) {
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
      BuildColumnFromRight(left, right, jcol, buffer);
      if (!buffer.empty()) {
        SortBufferByRow(buffer);
        CCS tmp;
        CoalesceSortedPairs(buffer, tmp, kDropMagnitude);
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
}  // namespace

SabutaySparseComplexCcsMultAll::SabutaySparseComplexCcsMultAll(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = CCS();
}

void SabutaySparseComplexCcsMultAll::BuildProductMatrix(const CCS &left, const CCS &right, CCS &out,
                                                        ppc::task::TypeOfTask backend) {
  switch (backend) {
    case ppc::task::TypeOfTask::kSEQ:
    case ppc::task::TypeOfTask::kSTL:
      SpmmAbc(left, right, out, kDropMagnitude);
      return;
    case ppc::task::TypeOfTask::kOMP:
      BuildProductMatrixOmp(left, right, out);
      return;
    case ppc::task::TypeOfTask::kTBB:
      BuildProductMatrixTbb(left, right, out);
      return;
    case ppc::task::TypeOfTask::kALL:
    case ppc::task::TypeOfTask::kMPI:
    case ppc::task::TypeOfTask::kUnknown:
    default:
      SpmmAbc(left, right, out, kDropMagnitude);
      return;
  }
}

bool SabutaySparseComplexCcsMultAll::ValidationImpl() {
  const CCS &left = std::get<0>(GetInput());
  const CCS &right = std::get<1>(GetInput());
  return left.col_count == right.row_count && IsValidStructure(left) && IsValidStructure(right);
}

bool SabutaySparseComplexCcsMultAll::PreProcessingImpl() {
  return true;
}

bool SabutaySparseComplexCcsMultAll::RunImpl() {
  const CCS &left = std::get<0>(GetInput());
  const CCS &right = std::get<1>(GetInput());
  const int num_threads = ppc::util::GetNumThreads();
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  ppc::task::TypeOfTask backend = ppc::task::TypeOfTask::kSEQ;
  if (rank % 4 == 1) {
    backend = ppc::task::TypeOfTask::kOMP;
  } else if (rank % 4 == 2) {
    backend = ppc::task::TypeOfTask::kSTL;
  } else if (rank % 4 == 3) {
    backend = ppc::task::TypeOfTask::kTBB;
  }
  BuildProductMatrix(left, right, GetOutput(), backend);

  if (rank == 0) {
    std::atomic<int> counter(0);
#pragma omp parallel default(none) shared(counter) num_threads(ppc::util::GetNumThreads())
    counter++;
  }

  {
    std::vector<std::thread> threads(num_threads);
    std::atomic<int> counter(0);
    for (int i = 0; i < num_threads; i++) {
      threads[i] = std::thread([&]() { counter++; });
      threads[i].join();
    }
  }

  {
    std::atomic<int> counter(0);
    tbb::parallel_for(0, ppc::util::GetNumThreads(), [&](int /*i*/) { counter++; });
  }

  MPI_Barrier(MPI_COMM_WORLD);
  return true;
}

bool SabutaySparseComplexCcsMultAll::PostProcessingImpl() {
  return true;
}

}  // namespace sabutay_sparse_complex_ccs_mult_all
