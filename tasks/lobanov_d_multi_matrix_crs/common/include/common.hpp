#pragma once

#include <cstddef>
#include <ostream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "task/include/task.hpp"

namespace lobanov_d_multi_matrix_crs {

struct CompressedRowMatrix {
  int row_count = 0;
  int column_count = 0;
  int non_zero_count = 0;
  std::vector<double> value_data;
  std::vector<int> column_index_data;
  std::vector<int> row_pointer_data;

  CompressedRowMatrix() = default;

  CompressedRowMatrix(int r, int c, int nz) : row_count(r), column_count(c), non_zero_count(nz) {
    if (nz > 0) {
      value_data.reserve(static_cast<std::size_t>(nz));
      column_index_data.reserve(static_cast<std::size_t>(nz));
    }
    row_pointer_data.reserve(static_cast<std::size_t>(r) + 1);
  }

  CompressedRowMatrix(const CompressedRowMatrix &other) = default;
  CompressedRowMatrix &operator=(const CompressedRowMatrix &other) = default;

  ~CompressedRowMatrix() = default;

  void ZeroInitialize() {
    row_count = 0;
    column_count = 0;
    non_zero_count = 0;
    value_data.clear();
    column_index_data.clear();
    row_pointer_data.clear();
  }

  [[nodiscard]] bool IsValid() const {
    if (row_count < 0 || column_count < 0 || non_zero_count < 0) {
      return false;
    }
    if (non_zero_count > 0) {
      if (value_data.size() != static_cast<std::size_t>(non_zero_count)) {
        return false;
      }
      if (column_index_data.size() != static_cast<std::size_t>(non_zero_count)) {
        return false;
      }
    }
    if (row_pointer_data.size() != static_cast<std::size_t>(row_count) + 1U) {
      return false;
    }
    if (!row_pointer_data.empty() && row_pointer_data[0] != 0) {
      return false;
    }
    return true;
  }
};

using InType = std::pair<CompressedRowMatrix, CompressedRowMatrix>;
using OutType = CompressedRowMatrix;
using TestType = std::tuple<std::string, CompressedRowMatrix, CompressedRowMatrix, CompressedRowMatrix>;
using BaseTask = ppc::task::Task<InType, OutType>;

CompressedRowMatrix CreateRandomCompressedRowMatrix(int row_count, int column_count, double density_factor,
                                                    int seed);

inline std::ostream &operator<<(std::ostream &os, const CompressedRowMatrix &matrix) {
  os << "CompressedRowMatrix{"
     << "rows=" << matrix.row_count << ", cols=" << matrix.column_count << ", nnz=" << matrix.non_zero_count << "}";
  return os;
}
}  // namespace lobanov_d_multi_matrix_crs