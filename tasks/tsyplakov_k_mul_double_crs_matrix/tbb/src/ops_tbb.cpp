#include "tsyplakov_k_mul_double_crs_matrix/tbb/include/ops_tbb.hpp"

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

#include <cmath>
#include <unordered_map>
#include <vector>

#include "tsyplakov_k_mul_double_crs_matrix/common/include/common.hpp"

namespace tsyplakov_k_mul_double_crs_matrix {

TsyplakovKTestTaskTBB::TsyplakovKTestTaskTBB(const InType& in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool TsyplakovKTestTaskTBB::ValidationImpl() {
  const auto& input = GetInput();
  return input.a.cols == input.b.rows;
}

bool TsyplakovKTestTaskTBB::PreProcessingImpl() { return true; }

bool TsyplakovKTestTaskTBB::RunImpl() {
  const auto& input = GetInput();
  const auto& A = input.a;
  const auto& B = input.b;

  const int rows = A.rows;

  std::vector<std::vector<double>> row_values(rows);
  std::vector<std::vector<int>> row_cols(rows);

  tbb::parallel_for(
      tbb::blocked_range<int>(0, rows),
      [&](const tbb::blocked_range<int>& r) {
        for (int i = r.begin(); i < r.end(); ++i) {
          std::unordered_map<int, double> acc;

          for (int a = A.row_ptr[i]; a < A.row_ptr[i + 1]; ++a) {
            const int k = A.col_index[a];
            const double valA = A.values[a];

            for (int b = B.row_ptr[k]; b < B.row_ptr[k + 1]; ++b) {
              const int j = B.col_index[b];
              acc[j] += valA * B.values[b];
            }
          }

          row_values[i].reserve(acc.size());
          row_cols[i].reserve(acc.size());

          for (const auto& [col, val] : acc) {
            if (std::fabs(val) > 1e-12) {
              row_cols[i].push_back(col);
              row_values[i].push_back(val);
            }
          }
        }
      });

  SparseMatrixCRS C(A.rows, B.cols);

  for (int i = 0; i < C.rows; ++i) {
    C.row_ptr[i + 1] =
        C.row_ptr[i] + static_cast<int>(row_values[i].size());
  }

  const int nnz = C.row_ptr[C.rows];

  C.values.reserve(nnz);
  C.col_index.reserve(nnz);

  for (int i = 0; i < C.rows; ++i) {
    C.values.insert(C.values.end(),
                    row_values[i].begin(),
                    row_values[i].end());

    C.col_index.insert(C.col_index.end(),
                       row_cols[i].begin(),
                       row_cols[i].end());
  }

  GetOutput() = std::move(C);

  return true;
}

bool TsyplakovKTestTaskTBB::PostProcessingImpl() { return true; }

}  // namespace tsyplakov_k_mul_double_crs_matrix