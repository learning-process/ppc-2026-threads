#pragma once

#include "akhmetov_daniil_strassen_dense_double_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace akhmetov_daniil_strassen_dense_double_seq {

class AkhmetovDStrassenDenseDoubleSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit AkhmetovDStrassenDenseDoubleSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace akhmetov_daniil_strassen_dense_double_seq
