#pragma once

#include "smetanin_v_hoare_even_odd_batchelor/common/include/common.hpp"
#include "task/include/task.hpp"

namespace smetanin_v_hoare_even_odd_batchelor {

class SmetaninVHoarSortSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit SmetaninVHoarSortSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace smetanin_v_hoare_even_odd_batchelor
