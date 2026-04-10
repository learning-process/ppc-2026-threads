#pragma once

#include <cstddef>
#include <vector>

#include "task/include/task.hpp"
#include "zyazeva_s_matrix_mult_cannon_alg/common/include/common.hpp"

#ifdef _OPENMP
#  include <omp.h>
#endif

namespace zyazeva_s_matrix_mult_cannon_alg {

class ZyazevaSMatrixMultCannonAlgOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit ZyazevaSMatrixMultCannonAlgOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static bool IsPerfectSquare(int x);
  static void MultiplyBlocks(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c,
                             int block_size);
  static void RegularMultiplication(const std::vector<double> &m1, const std::vector<double> &m2,
                                    std::vector<double> &res, int sz);
  static void InitializeBlocks(const std::vector<double> &m1, const std::vector<double> &m2,
                               std::vector<std::vector<double>> &blocks_a, std::vector<std::vector<double>> &blocks_b,
                               int grid_size, int block_size, size_t grid_size_t, size_t block_size_t, size_t sz_t);
  static void AlignBlocks(const std::vector<std::vector<double>> &blocks_a,
                          const std::vector<std::vector<double>> &blocks_b, std::vector<std::vector<double>> &aligned_a,
                          std::vector<std::vector<double>> &aligned_b, int grid_size, size_t grid_size_t);
  static void CannonStep(std::vector<std::vector<double>> &aligned_a, std::vector<std::vector<double>> &aligned_b,
                         std::vector<std::vector<double>> &blocks_c, int grid_size, int block_size, size_t grid_size_t,
                         int step);
  static void AssembleResult(const std::vector<std::vector<double>> &blocks_c, std::vector<double> &res_m,
                             int grid_size, int block_size, size_t sz_t, size_t grid_size_t, size_t block_size_t);
};

}  // namespace zyazeva_s_matrix_mult_cannon_alg
