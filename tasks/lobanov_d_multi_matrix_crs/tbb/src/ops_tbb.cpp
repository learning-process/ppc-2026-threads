#include "lobanov_d_multi_matrix_crs/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "lobanov_d_multi_matrix_crs/common/include/common.hpp"

namespace lobanov_d_multi_matrix_crs {

namespace {

void ProcessRow(const CompressedRowMatrix &a, const CompressedRowMatrix &b, int row_idx,
                std::vector<std::pair<int, double>> &row_pairs) {
  row_pairs.clear();

  for (int pos_a = a.row_pointer_data[row_idx]; pos_a < a.row_pointer_data[row_idx + 1]; ++pos_a) {
    int col_a = a.column_index_data[pos_a];
    double val_a = a.value_data[pos_a];

    for (int pos_b = b.row_pointer_data[col_a]; pos_b < b.row_pointer_data[col_a + 1]; ++pos_b) {
      int col_b = b.column_index_data[pos_b];
      double val_b = b.value_data[pos_b];
      row_pairs.emplace_back(col_b, val_a * val_b);
    }
  }
}

void MergeRowPairs(std::vector<std::pair<int, double>> &row_pairs) {
  if (row_pairs.empty()) {
    return;
  }

  std::ranges::sort(row_pairs, [](const auto &lhs, const auto &rhs) { return lhs.first < rhs.first; });

  std::vector<std::pair<int, double>> merged;
  double sum = row_pairs[0].second;
  int prev_col = row_pairs[0].first;

  for (std::size_t k = 1; k < row_pairs.size(); ++k) {
    if (row_pairs[k].first == prev_col) {
      sum += row_pairs[k].second;
    } else {
      if (std::abs(sum) > 1e-15) {
        merged.emplace_back(prev_col, sum);
      }
      prev_col = row_pairs[k].first;
      sum = row_pairs[k].second;
    }
  }
  if (std::abs(sum) > 1e-15) {
    merged.emplace_back(prev_col, sum);
  }

  row_pairs = std::move(merged);
}

}  // namespace

LobanovMultyMatrixTBB::LobanovMultyMatrixTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool LobanovMultyMatrixTBB::ValidationImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());
  return a.column_count == b.row_count && a.row_count > 0 && b.column_count > 0;
}

bool LobanovMultyMatrixTBB::PreProcessingImpl() {
  GetOutput() = {};
  return true;
}

bool LobanovMultyMatrixTBB::RunImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());
  auto &c = GetOutput();

  c.row_count = a.row_count;
  c.column_count = b.column_count;
  c.row_pointer_data.assign(static_cast<std::size_t>(c.row_count) + 1U, 0);
  c.non_zero_count = 0;
  c.value_data.clear();
  c.column_index_data.clear();

  if (a.row_count == 0 || b.column_count == 0 || a.non_zero_count == 0 || b.non_zero_count == 0) {
    return true;
  }

  std::vector<std::vector<std::pair<int, double>>> local_entries(static_cast<std::size_t>(a.row_count));

  tbb::parallel_for(tbb::blocked_range<int>(0, a.row_count), [&](const tbb::blocked_range<int> &range) {
    for (int i = range.begin(); i != range.end(); ++i) {
      std::vector<std::pair<int, double>> row_pairs;
      ProcessRow(a, b, i, row_pairs);
      MergeRowPairs(row_pairs);
      local_entries[static_cast<std::size_t>(i)] = std::move(row_pairs);
    }
  });

  for (int i = 0; i < c.row_count; ++i) {
    c.row_pointer_data[i + 1] =
        c.row_pointer_data[i] + static_cast<int>(local_entries[static_cast<std::size_t>(i)].size());
  }

  const int total_nnz = c.row_pointer_data[c.row_count];
  c.value_data.reserve(static_cast<std::size_t>(total_nnz));
  c.column_index_data.reserve(static_cast<std::size_t>(total_nnz));

  for (int i = 0; i < c.row_count; ++i) {
    for (const auto &[col, val] : local_entries[static_cast<std::size_t>(i)]) {
      c.value_data.push_back(val);
      c.column_index_data.push_back(col);
    }
  }
  c.non_zero_count = total_nnz;

  return true;
}

bool LobanovMultyMatrixTBB::PostProcessingImpl() {
  return true;
}

}  // namespace lobanov_d_multi_matrix_crs
