#pragma once

#include "alekseev_a_mult_matrix_crs/common/include/common.hpp"
#include "task/include/task.hpp"

namespace alekseev_a_mult_matrix_crs {

class AlekseevAMultMatrixCRSSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }

  explicit AlekseevAMultMatrixCRSSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace alekseev_a_mult_matrix_crs
