#include "sabutay_sparse_complex_ccs_mult_all/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <utility>

#include "sabutay_sparse_complex_ccs_mult_all/all/include/ops_all.hpp"
#include "sabutay_sparse_complex_ccs_mult_all/common/include/common.hpp"
#include "task/include/task.hpp"

namespace sabutay_sparse_complex_ccs_mult_all {
namespace {
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
}  // namespace

SabutaySparseComplexCcsMultFixSEQ::SabutaySparseComplexCcsMultFixSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = CCS();
}

void SabutaySparseComplexCcsMultFixSEQ::BuildProductMatrix(const CCS &left, const CCS &right, CCS &out) {
  SabutaySparseComplexCcsMultAll::BuildProductMatrix(left, right, out, ppc::task::TypeOfTask::kSEQ);
}

bool SabutaySparseComplexCcsMultFixSEQ::ValidationImpl() {
  const CCS &left = std::get<0>(GetInput());
  const CCS &right = std::get<1>(GetInput());
  return left.col_count == right.row_count && IsValidStructure(left) && IsValidStructure(right);
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

}  // namespace sabutay_sparse_complex_ccs_mult_all
