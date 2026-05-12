#pragma once

#include <vector>

#include "maslova_u_mult_matr_crs/common/include/common.hpp"
#include "task/include/task.hpp"

namespace maslova_u_mult_matr_crs {

class MaslovaUMultMatrALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit MaslovaUMultMatrALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
  static void BroadcastCRSMatrix(CRSMatrix &matrix, int root, int rows, int cols);

  static void ComputeLocalPart(const CRSMatrix &a, const CRSMatrix &b, int start_row, int local_rows,
                               std::vector<int> &local_nnz, std::vector<double> &flat_values,
                               std::vector<int> &flat_cols);

  static void GatherResults(int rank, int size, int a_rows, int b_cols, int local_rows, CRSMatrix &c,
                            const std::vector<int> &local_nnz, const std::vector<double> &flat_values,
                            const std::vector<int> &flat_cols);
  static void SortVector(std::vector<int> &vec);
};

}  // namespace maslova_u_mult_matr_crs
