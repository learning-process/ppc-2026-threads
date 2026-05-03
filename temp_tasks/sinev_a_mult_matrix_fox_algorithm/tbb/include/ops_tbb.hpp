#pragma once

#include <cstddef>
#include <vector>

#include "sinev_a_mult_matrix_fox_algorithm/common/include/common.hpp"
#include "task/include/task.hpp"

namespace sinev_a_mult_matrix_fox_algorithm {

class SinevAMultMatrixFoxAlgorithmTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit SinevAMultMatrixFoxAlgorithmTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void SimpleMultiply(size_t n, const std::vector<double> &a, const std::vector<double> &b,
                             std::vector<double> &c);
  static void DecomposeToBlocks(const std::vector<double> &src, std::vector<double> &dst, size_t n, size_t bs, int q);
  static void AssembleFromBlocks(const std::vector<double> &src, std::vector<double> &dst, size_t n, size_t bs, int q);
  static void FoxStep(const std::vector<double> &blocks_a, const std::vector<double> &blocks_b,
                      std::vector<double> &blocks_c, size_t bs, int q, int step);
  static void MultiplyBlocks(const std::vector<double> &blocks_a, const std::vector<double> &blocks_b,
                             std::vector<double> &blocks_c, size_t bs, size_t a_off, size_t b_off, size_t c_off);
  static size_t ChooseBlockSize(size_t n);
};

}  // namespace sinev_a_mult_matrix_fox_algorithm
