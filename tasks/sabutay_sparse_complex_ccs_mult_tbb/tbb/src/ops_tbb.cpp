#include "sabutay_sparse_complex_ccs_mult_tbb/tbb/include/ops_tbb.hpp"

#include <algorithm>
#include <cstddef>
#include <vector>

#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/enumerable_thread_specific.h"
#include "oneapi/tbb/parallel_for.h"
#include "sabutay_sparse_complex_ccs_mult_tbb/common/include/common.hpp"

namespace sabutay_sparse_complex_ccs_mult_tbb {

namespace {

struct ColumnData {
  std::vector<int> row_ind;
  std::vector<Complex> values;
};

struct ThreadWorkspace {
  explicit ThreadWorkspace(int rows)
      : accumulator(static_cast<std::size_t>(rows), Complex{}), marker(static_cast<std::size_t>(rows), -1) {}

  std::vector<Complex> accumulator;
  std::vector<int> marker;
  std::vector<int> touched_rows;
};

}  // namespace

SabutayASparseComplexCcsMultTBB::SabutayASparseComplexCcsMultTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = SparseMatrixCCS{};
}

SparseMatrixCCS SabutayASparseComplexCcsMultTBB::MultiplyTbb(const SparseMatrixCCS &lhs, const SparseMatrixCCS &rhs) {
  if (!IsValidCcs(lhs) || !IsValidCcs(rhs) || lhs.cols != rhs.rows) {
    return {};
  }

  SparseMatrixCCS result = MakeZeroMatrix(lhs.rows, rhs.cols);
  std::vector<ColumnData> columns(static_cast<std::size_t>(rhs.cols));

  oneapi::tbb::enumerable_thread_specific<ThreadWorkspace> workspaces([&lhs]() { return ThreadWorkspace(lhs.rows); });

  oneapi::tbb::parallel_for(oneapi::tbb::blocked_range<int>(0, rhs.cols),
                            [&](const oneapi::tbb::blocked_range<int> &range) {
    ThreadWorkspace &workspace = workspaces.local();

    for (int col = range.begin(); col < range.end(); ++col) {
      workspace.touched_rows.clear();

      for (int rhs_idx = rhs.col_ptr[static_cast<std::size_t>(col)];
           rhs_idx < rhs.col_ptr[static_cast<std::size_t>(col) + 1U]; ++rhs_idx) {
        const int lhs_col = rhs.row_ind[static_cast<std::size_t>(rhs_idx)];
        const Complex rhs_value = rhs.values[static_cast<std::size_t>(rhs_idx)];

        for (int lhs_idx = lhs.col_ptr[static_cast<std::size_t>(lhs_col)];
             lhs_idx < lhs.col_ptr[static_cast<std::size_t>(lhs_col) + 1U]; ++lhs_idx) {
          const int row = lhs.row_ind[static_cast<std::size_t>(lhs_idx)];

          if (workspace.marker[static_cast<std::size_t>(row)] != col) {
            workspace.marker[static_cast<std::size_t>(row)] = col;
            workspace.accumulator[static_cast<std::size_t>(row)] = Complex{};
            workspace.touched_rows.push_back(row);
          }

          workspace.accumulator[static_cast<std::size_t>(row)] +=
              lhs.values[static_cast<std::size_t>(lhs_idx)] * rhs_value;
        }
      }

      std::sort(workspace.touched_rows.begin(), workspace.touched_rows.end());

      ColumnData &column_data = columns[static_cast<std::size_t>(col)];
      column_data.row_ind.clear();
      column_data.values.clear();
      column_data.row_ind.reserve(workspace.touched_rows.size());
      column_data.values.reserve(workspace.touched_rows.size());

      for (const int row : workspace.touched_rows) {
        const Complex value = workspace.accumulator[static_cast<std::size_t>(row)];
        if (!IsNearZero(value)) {
          column_data.row_ind.push_back(row);
          column_data.values.push_back(value);
        }
      }
    }
  });

  int nnz = 0;
  for (int col = 0; col < rhs.cols; ++col) {
    result.col_ptr[static_cast<std::size_t>(col)] = nnz;
    nnz += static_cast<int>(columns[static_cast<std::size_t>(col)].values.size());
  }
  result.col_ptr[static_cast<std::size_t>(rhs.cols)] = nnz;

  result.row_ind.reserve(static_cast<std::size_t>(nnz));
  result.values.reserve(static_cast<std::size_t>(nnz));

  for (int col = 0; col < rhs.cols; ++col) {
    const ColumnData &column_data = columns[static_cast<std::size_t>(col)];
    result.row_ind.insert(result.row_ind.end(), column_data.row_ind.begin(), column_data.row_ind.end());
    result.values.insert(result.values.end(), column_data.values.begin(), column_data.values.end());
  }

  return result;
}

bool SabutayASparseComplexCcsMultTBB::ValidationImpl() {
  const auto &[lhs, rhs] = GetInput();
  is_input_valid_ = IsValidCcs(lhs) && IsValidCcs(rhs) && lhs.cols == rhs.rows;
  return is_input_valid_;
}

bool SabutayASparseComplexCcsMultTBB::PreProcessingImpl() {
  if (!is_input_valid_) {
    lhs_ = SparseMatrixCCS{};
    rhs_ = SparseMatrixCCS{};
    result_ = SparseMatrixCCS{};
    GetOutput() = result_;
    return true;
  }

  const auto &[lhs, rhs] = GetInput();
  lhs_ = NormalizeCcs(lhs);
  rhs_ = NormalizeCcs(rhs);
  result_ = MakeZeroMatrix(lhs_.rows, rhs_.cols);
  GetOutput() = result_;
  return true;
}

bool SabutayASparseComplexCcsMultTBB::RunImpl() {
  if (!is_input_valid_) {
    return true;
  }

  result_ = MultiplyTbb(lhs_, rhs_);
  return IsValidCcs(result_);
}

bool SabutayASparseComplexCcsMultTBB::PostProcessingImpl() {
  GetOutput() = result_;
  return true;
}

}  // namespace sabutay_sparse_complex_ccs_mult_tbb
