#pragma once

#include <string>

#include "task/include/task.hpp"
#include "yakimov_i_mult_of_dense_matrices_fox_algorithm_seq/common/include/common.hpp"

namespace yakimov_i_mult_of_dense_matrices_fox_algorithm_seq {

class YakimovIMultOfDenseMatricesFoxAlgorithmOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit YakimovIMultOfDenseMatricesFoxAlgorithmOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  DenseMatrix matrix_a_;
  DenseMatrix matrix_b_;
  DenseMatrix result_matrix_;
  std::string matrix_a_filename_;
  std::string matrix_b_filename_;
  int block_size_ = 0;
};

}  // namespace yakimov_i_mult_of_dense_matrices_fox_algorithm_seq
