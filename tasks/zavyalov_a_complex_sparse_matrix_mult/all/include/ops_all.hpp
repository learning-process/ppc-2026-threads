#pragma once

#include <cstddef>
#include <map>
#include <utility>

#include "task/include/task.hpp"
#include "zavyalov_a_complex_sparse_matrix_mult/common/include/common.hpp"

namespace zavyalov_a_compl_sparse_matr_mult {

class ZavyalovAComplSparseMatrMultALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit ZavyalovAComplSparseMatrMultALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
  static std::map<std::pair<size_t, size_t>, Complex> ComputeLocalChunk(const SparseMatrix &matr_a,
                                                                        const SparseMatrix &matr_b, size_t start,
                                                                        size_t end);
  SparseMatrix MultiplicateWithMPI(const SparseMatrix &matr_a, const SparseMatrix &matr_b);
};

}  // namespace zavyalov_a_compl_sparse_matr_mult
