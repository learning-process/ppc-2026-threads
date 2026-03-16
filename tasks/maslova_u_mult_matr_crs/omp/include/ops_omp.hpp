#pragma once

#include <vector>

#include "maslova_u_mult_matr_crs/common/include/common.hpp"
#include "task/include/task.hpp"

namespace maslova_u_mult_matr_crs {

class MaslovaUMultMatrOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit MaslovaUMultMatrOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
  static int GetRowNNZ(int i, const CRSMatrix &a, const CRSMatrix &b, std::vector<int> &marker);
  static void FillRowValues(int i, const CRSMatrix &a, const CRSMatrix &b, CRSMatrix &c, std::vector<double> &acc,
                            std::vector<int> &marker, std::vector<int> &used);
};

}  // namespace maslova_u_mult_matr_crs
