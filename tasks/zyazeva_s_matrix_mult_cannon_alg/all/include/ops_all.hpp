#pragma once

#include <cstddef>
#include <vector>

#include "task/include/task.hpp"
#include "zyazeva_s_matrix_mult_cannon_alg/common/include/common.hpp"

#ifdef _OPENMP
#  include <omp.h>
#endif

namespace zyazeva_s_matrix_mult_cannon_alg {

class ZyazevaSMatrixMultCannonAlgALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit ZyazevaSMatrixMultCannonAlgALL(const InType &in);

 private:
  int rank_{0};
  int mpi_size_{1};

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

  // extracted to reduce RunImpl cognitive complexity
  void DistributeBlocks(const std::vector<double> &m1, const std::vector<double> &m2, int grid, int block_size,
                        int block_elems, int sz, std::vector<double> &local_a, std::vector<double> &local_b) const;
  void CollectResult(const std::vector<double> &local_c, std::vector<double> &result, int grid, int block_size,
                     int block_elems, int sz) const;
};

}  // namespace zyazeva_s_matrix_mult_cannon_alg
