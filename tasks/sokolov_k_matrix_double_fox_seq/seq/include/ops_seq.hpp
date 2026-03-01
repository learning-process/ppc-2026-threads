#pragma once

#include <vector>

#include "sokolov_k_matrix_double_fox_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace sokolov_k_matrix_double_fox_seq {

class SokolovKMatrixDoubleFoxSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit SokolovKMatrixDoubleFoxSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  int n_{};
  int block_size_{};
  int q_{};
  std::vector<double> blocks_a_;
  std::vector<double> blocks_b_;
  std::vector<double> blocks_c_;
};

}  // namespace sokolov_k_matrix_double_fox_seq
