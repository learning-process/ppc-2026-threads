#pragma once

#include <vector>

#include "maslova_u_mult_matr_crs/common/include/common.hpp"
#include "task/include/task.hpp"

namespace maslova_u_mult_matr_crs {

class MaslovaUMultMatrSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit MaslovaUMultMatrSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  void ProcessRow(int i, const CRSMatrix &a, const CRSMatrix &b, CRSMatrix &c);

  std::vector<double> temp_row_;
  std::vector<int> marker_;
  std::vector<int> used_cols_;
};

}  // namespace maslova_u_mult_matr_crs
