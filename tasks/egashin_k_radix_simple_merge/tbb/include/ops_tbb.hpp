#pragma once

#include "egashin_k_radix_simple_merge/common/include/common.hpp"
#include "task/include/task.hpp"

namespace egashin_k_radix_simple_merge {

class EgashinKRadixSimpleMergeTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }

  explicit EgashinKRadixSimpleMergeTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  OutType result_;
};

}  // namespace egashin_k_radix_simple_merge
