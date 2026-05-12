#include "posternak_a_crs_mul_complex_matrix/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

#include "posternak_a_crs_mul_complex_matrix/common/include/common.hpp"

namespace {

size_t ComputeRowNoZeroCount(const posternak_a_crs_mul_complex_matrix::CRSMatrix &a,
                             const posternak_a_crs_mul_complex_matrix::CRSMatrix &b, int row, double threshold) {
  std::unordered_map<int, std::complex<double>> row_sum;

  for (int idx_a = a.index_row[row]; idx_a < a.index_row[row + 1]; ++idx_a) {
    int col_a = a.index_col[idx_a];
    auto val_a = a.values[idx_a];

    for (int idx_b = b.index_row[col_a]; idx_b < b.index_row[col_a + 1]; ++idx_b) {
      int col_b = b.index_col[idx_b];
      auto val_b = b.values[idx_b];
      row_sum[col_b] += val_a * val_b;
    }
  }

  size_t local = 0;
  for (const auto &[col, val] : row_sum) {
    if (std::abs(val) > threshold) {
      ++local;
    }
  }
  return local;
}

void BuildResultStructure(posternak_a_crs_mul_complex_matrix::CRSMatrix &res, std::vector<size_t> &row_prefix) {
  for (int i = 1; i < res.rows; ++i) {
    row_prefix[i] += row_prefix[i - 1];
  }

  const size_t total = row_prefix.empty() ? 0 : row_prefix.back();
  res.values.resize(total);
  res.index_col.resize(total);
  res.index_row.resize(res.rows + 1);

  for (int i = 0; i <= res.rows; ++i) {
    res.index_row[i] = (i == 0 ? 0 : static_cast<int>(row_prefix[i - 1]));
  }
}

void ComputeAndWriteRow(const posternak_a_crs_mul_complex_matrix::CRSMatrix &a,
                        const posternak_a_crs_mul_complex_matrix::CRSMatrix &b,
                        posternak_a_crs_mul_complex_matrix::CRSMatrix &res, int row, double threshold) {
  std::unordered_map<int, std::complex<double>> row_sum;

  for (int idx_a = a.index_row[row]; idx_a < a.index_row[row + 1]; ++idx_a) {
    int col_a = a.index_col[idx_a];
    auto val_a = a.values[idx_a];

    for (int idx_b = b.index_row[col_a]; idx_b < b.index_row[col_a + 1]; ++idx_b) {
      int col_b = b.index_col[idx_b];
      auto val_b = b.values[idx_b];
      row_sum[col_b] += val_a * val_b;
    }
  }

  std::vector<std::pair<int, std::complex<double>>> sorted(row_sum.begin(), row_sum.end());

  std::ranges::sort(sorted, [](const auto &p1, const auto &p2) { return p1.first < p2.first; });

  size_t pos = res.index_row[row];
  for (const auto &[col_idx, value] : sorted) {
    if (std::abs(value) > threshold) {
      res.values[pos] = value;
      res.index_col[pos] = col_idx;
      ++pos;
    }
  }
}

bool HandleEmptyInput(posternak_a_crs_mul_complex_matrix::CRSMatrix &res) {
  res.values.clear();
  res.index_col.clear();
  res.index_row.assign(res.rows + 1, 0);
  return true;
}

std::vector<size_t> CountNonZeroElementsParallel(const posternak_a_crs_mul_complex_matrix::CRSMatrix &a,
                                                 const posternak_a_crs_mul_complex_matrix::CRSMatrix &b, int total_rows,
                                                 double threshold) {
  std::vector<size_t> no_zero_rows(total_rows);
  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 1;
  }

  const unsigned int chunk_size = std::max(1U, (static_cast<unsigned int>(total_rows) + num_threads - 1) / num_threads);

  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  for (unsigned int thr = 0; thr < num_threads; ++thr) {
    const int start_row = static_cast<int>(thr * chunk_size);
    const int end_row = static_cast<int>(
        std::min(static_cast<unsigned int>(start_row) + chunk_size, static_cast<unsigned int>(total_rows)));
    if (start_row >= total_rows) {
      break;
    }

    threads.emplace_back([&, start_row, end_row]() {
      for (int row = start_row; row < end_row; ++row) {
        no_zero_rows[row] = ComputeRowNoZeroCount(a, b, row, threshold);
      }
    });
  }

  for (auto &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }

  return no_zero_rows;
}

void ComputeResultValuesParallel(const posternak_a_crs_mul_complex_matrix::CRSMatrix &a,
                                 const posternak_a_crs_mul_complex_matrix::CRSMatrix &b,
                                 posternak_a_crs_mul_complex_matrix::CRSMatrix &res, int total_rows, double threshold) {
  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 1;
  }

  const unsigned int chunk_size = std::max(1U, (static_cast<unsigned int>(total_rows) + num_threads - 1) / num_threads);

  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  for (unsigned int thr = 0; thr < num_threads; ++thr) {
    const int start_row = static_cast<int>(thr * chunk_size);
    const int end_row = static_cast<int>(
        std::min(static_cast<unsigned int>(start_row) + chunk_size, static_cast<unsigned int>(total_rows)));
    if (start_row >= total_rows) {
      break;
    }

    threads.emplace_back([&, start_row, end_row]() {
      for (int row = start_row; row < end_row; ++row) {
        ComputeAndWriteRow(a, b, res, row, threshold);
      }
    });
  }

  for (auto &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }
}

}  // namespace

namespace posternak_a_crs_mul_complex_matrix {

PosternakACRSMulComplexMatrixSTL::PosternakACRSMulComplexMatrixSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = CRSMatrix{};
}

bool PosternakACRSMulComplexMatrixSTL::ValidationImpl() {
  const auto &input = GetInput();
  const auto &a = input.first;
  const auto &b = input.second;
  return a.IsValid() && b.IsValid() && a.cols == b.rows;
}

bool PosternakACRSMulComplexMatrixSTL::PreProcessingImpl() {
  const auto &input = GetInput();
  const auto &a = input.first;
  const auto &b = input.second;
  auto &res = GetOutput();

  res.rows = a.rows;
  res.cols = b.cols;
  return true;
}

bool PosternakACRSMulComplexMatrixSTL::RunImpl() {
  const auto &input = GetInput();
  const auto &a = input.first;
  const auto &b = input.second;
  auto &res = GetOutput();

  if (a.values.empty() || b.values.empty()) {
    return HandleEmptyInput(res);
  }

  constexpr double kThreshold = 1e-12;
  res.rows = a.rows;
  res.cols = b.cols;

  std::vector<size_t> no_zero_rows = CountNonZeroElementsParallel(a, b, res.rows, kThreshold);

  BuildResultStructure(res, no_zero_rows);

  ComputeResultValuesParallel(a, b, res, res.rows, kThreshold);

  return res.IsValid();
}

bool PosternakACRSMulComplexMatrixSTL::PostProcessingImpl() {
  return GetOutput().IsValid();
}

}  // namespace posternak_a_crs_mul_complex_matrix
