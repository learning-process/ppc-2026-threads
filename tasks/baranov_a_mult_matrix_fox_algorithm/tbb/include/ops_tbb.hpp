#pragma once

#include <cstddef>

#include "baranov_a_mult_matrix_fox_algorithm/common/include/common.hpp"
#include "task/include/task.hpp"

namespace baranov_a_mult_matrix_fox_algorithm_tbb {

class BaranovAMultMatrixFoxAlgorithmTBB : public baranov_a_mult_matrix_fox_algorithm::BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit BaranovAMultMatrixFoxAlgorithmTBB(const baranov_a_mult_matrix_fox_algorithm::InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  void FoxBlockMultiplication(size_t n, size_t block_size);
  void StandardMultiplication(size_t n);
};

}  // namespace baranov_a_mult_matrix_fox_algorithm_tbb
