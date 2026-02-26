#pragma once

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

  std::vector<double> temp_row;
  std::vector<int> marker;
  std::vector<int> used_cols;
};

}  // namespace maslova_u_mult_matr_crs
