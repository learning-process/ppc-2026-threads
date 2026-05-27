#pragma once

#include <vector>

#include "maslova_u_mult_matr_crs/common/include/common.hpp"
#include "task/include/task.hpp"

namespace maslova_u_mult_matr_crs {

class MaslovaUMultMatrSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit MaslovaUMultMatrSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
  void ComputeRowPtrSTL(int rows_a, int cols_b, int num_threads);
  void ComputeValuesSTL(int rows_a, int cols_b, int num_threads);
  static int GetRowNNZ(int i, const CRSMatrix &a, const CRSMatrix &b, std::vector<int> &marker);
  static void FillRowValues(int i, const CRSMatrix &a, const CRSMatrix &b, CRSMatrix &c, std::vector<double> &acc,
                            std::vector<int> &marker, std::vector<int> &used);
};

}  // namespace maslova_u_mult_matr_crs
