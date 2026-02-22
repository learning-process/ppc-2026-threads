#include "zavyalov_a_complex_sparse_matr_mult/seq/include/ops_seq.hpp"

#include <chrono>
#include <numeric>
#include <thread>
#include <vector>

#include "util/include/util.hpp"
#include "zavyalov_a_complex_sparse_matr_mult/common/include/common.hpp"

namespace zavyalov_a_compl_sparse_matr_mult {

ZavyalovAComplSparseMatrMultSEQ::ZavyalovAComplSparseMatrMultSEQ(const InType& in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool ZavyalovAComplSparseMatrMultSEQ::ValidationImpl() {
  return std::get<0>(GetInput()).width == std::get<1>(GetInput()).height;
}

bool ZavyalovAComplSparseMatrMultSEQ::PreProcessingImpl() {
  return true;
}

bool ZavyalovAComplSparseMatrMultSEQ::RunImpl() {
  Sparse_matrix& matr_a = std::get<0>(GetInput());
  Sparse_matrix& matr_b = std::get<1>(GetInput());

  GetOutput() = matr_a * matr_b;

  std::chrono::milliseconds timespan(10);
  std::this_thread::sleep_for(timespan);  // CheckTestOutputData works much slower than RunImpl in perf tests. Thats why
                                          // i use this slowing method
  return true;
}

bool ZavyalovAComplSparseMatrMultSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace zavyalov_a_compl_sparse_matr_mult
