#pragma once

#include <tbb/tbb.h>

#include "task/include/task.hpp"
#include "zyazeva_s_matrix_mult_cannon_alg/common/include/common.hpp"

namespace zyazeva_s_matrix_mult_cannon_alg {

class ZyazevaSMatrixMultCannonAlgTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }

  explicit ZyazevaSMatrixMultCannonAlgTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static void MultiplyBlocks(const double *a, const double *b, double *c, int block_size);
};

}  // namespace zyazeva_s_matrix_mult_cannon_alg
