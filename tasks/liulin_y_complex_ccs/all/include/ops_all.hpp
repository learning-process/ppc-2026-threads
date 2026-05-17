#pragma once

#include <mpi.h>

#include <complex>
#include <vector>

#include "liulin_y_complex_ccs/common/include/common.hpp"
#include "task/include/task.hpp"

namespace liulin_y_complex_ccs {

class LiulinYComplexCcsAll : public ppc::task::Task<InType, OutType> {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit LiulinYComplexCcsAll(const InType &in);

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

 private:
  void BroadcastMatrix(CCSMatrix &mat, int root, MPI_Comm comm);
  void Transpose(const CCSMatrix &in, CCSMatrix &out);
};

}  // namespace liulin_y_complex_ccs
