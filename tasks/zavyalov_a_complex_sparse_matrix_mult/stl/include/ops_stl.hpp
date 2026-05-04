#pragma once

#include <cstddef>
#include <map>
#include <utility>
#include <vector>

#include "task/include/task.hpp"
#include "zavyalov_a_complex_sparse_matrix_mult/common/include/common.hpp"

namespace zavyalov_a_compl_sparse_matr_mult {

class ZavyalovAComplSparseMatrMultSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit ZavyalovAComplSparseMatrMultSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
  static SparseMatrix MultiplicateWithStl(const SparseMatrix &matr_a, const SparseMatrix &matr_b);
  static void Worker(int tid, std::size_t start, std::size_t end, const SparseMatrix &matr_a,
                     const SparseMatrix &matr_b,
                     std::vector<std::map<std::pair<std::size_t, std::size_t>, Complex>> &local_maps);
};

}  // namespace zavyalov_a_compl_sparse_matr_mult
