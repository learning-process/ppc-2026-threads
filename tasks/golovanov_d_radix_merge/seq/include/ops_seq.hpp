#pragma once

#include "golovanov_d_radix_merge/common/include/common.hpp"
#include "task/include/task.hpp"

namespace golovanov_d_radix_merge {

class GolovanovDRadixMergeSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit GolovanovDRadixMergeSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace golovanov_d_radix_merge
