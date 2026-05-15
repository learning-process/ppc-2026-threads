#pragma once

#include "ermakov_a_spar_mat_mult/common/include/common.hpp"
#include "task/include/task.hpp"

namespace ermakov_a_spar_mat_mult {

class ErmakovASparMatMultALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit ErmakovASparMatMultALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static bool ValidateMatrix(const MatrixCRS &m);

  MatrixCRS a_;
  MatrixCRS b_;
  MatrixCRS c_;
};

}  // namespace ermakov_a_spar_mat_mult
