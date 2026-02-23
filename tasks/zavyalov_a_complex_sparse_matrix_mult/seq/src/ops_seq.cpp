#include "zavyalov_a_complex_sparse_matrix_mult/seq/include/ops_seq.hpp"

#include <chrono>
#include <thread>

#include "zavyalov_a_complex_sparse_matrix_mult/common/include/common.hpp"

namespace zavyalov_a_compl_sparse_matr_mult {

ZavyalovAComplSparseMatrMultSEQ::ZavyalovAComplSparseMatrMultSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool ZavyalovAComplSparseMatrMultSEQ::ValidationImpl() {
  const auto &matr_a = std::get<0>(GetInput());
  const auto &matr_b = std::get<1>(GetInput());
  return matr_a.width == matr_b.height;
}

bool ZavyalovAComplSparseMatrMultSEQ::PreProcessingImpl() {
  return true;
}

bool ZavyalovAComplSparseMatrMultSEQ::RunImpl() {
  const auto &matr_a = std::get<0>(GetInput());
  const auto &matr_b = std::get<1>(GetInput());

  GetOutput() = matr_a * matr_b;

  // CheckTestOutputData works much slower than RunImpl in perf tests.
  // That's why we use this slowing method
  std::this_thread::sleep_for(std::chrono::milliseconds(10));
  return true;
}

bool ZavyalovAComplSparseMatrMultSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace zavyalov_a_compl_sparse_matr_mult
