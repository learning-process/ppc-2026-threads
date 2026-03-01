#include "luzan_e_double_sparse_matrix_mult_seq/seq/include/ops_seq.hpp"

#include <numeric>
#include <vector>

#include "luzan_e_double_sparse_matrix_mult_seq/common/include/common.hpp"
#include "util/include/util.hpp"

namespace luzan_e_double_sparse_matrix_mult_seq {

LuzanEDoubleSparseMatrixMultSeq::LuzanEDoubleSparseMatrixMultSeq(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  // GetOutput() = 0;
}

bool LuzanEDoubleSparseMatrixMultSeq::ValidationImpl() {
  const auto &A = std::get<0>(GetInput());
  const auto &B = std::get<1>(GetInput());
  return A.GetCols() == B.GetRows() && 
        A.GetCols() != 0 && A.GetRows() != 0 && B.GetCols() != 0;
}

bool LuzanEDoubleSparseMatrixMultSeq::PreProcessingImpl() {
  return true;
}

bool LuzanEDoubleSparseMatrixMultSeq::RunImpl() {
  const auto &A = std::get<0>(GetInput());
  const auto &B = std::get<1>(GetInput());

  GetOutput() = A * B;
  return true;
}

bool LuzanEDoubleSparseMatrixMultSeq::PostProcessingImpl() {
  return true;
}

}  // namespace luzan_e_double_sparse_matrix_mult_seq
