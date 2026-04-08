#include "luzan_e_double_sparse_matrix_mult_seq/seq/include/ops_seq.hpp"

#include "luzan_e_double_sparse_matrix_mult_seq/common/include/common.hpp"
// #include "util/include/util.hpp"

namespace luzan_e_double_sparse_matrix_mult_seq {

LuzanEDoubleSparseMatrixMultSeq::LuzanEDoubleSparseMatrixMultSeq(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  // GetOutput() = 0;
}

bool LuzanEDoubleSparseMatrixMultSeq::ValidationImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());
  return a.GetCols() == b.GetRows() && a.GetCols() != 0 && a.GetRows() != 0 && b.GetCols() != 0;
}

bool LuzanEDoubleSparseMatrixMultSeq::PreProcessingImpl() {
  return true;
}

bool LuzanEDoubleSparseMatrixMultSeq::RunImpl() {
  const auto &a = std::get<0>(GetInput());
  const auto &b = std::get<1>(GetInput());

  GetOutput() = a * b;
  return true;
}

bool LuzanEDoubleSparseMatrixMultSeq::PostProcessingImpl() {
  return true;
}

}  // namespace luzan_e_double_sparse_matrix_mult_seq
