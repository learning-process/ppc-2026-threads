#pragma once

#include "smetanin_d_hoare_even_odd_batchelor/common/include/common.hpp"
#include "task/include/task.hpp"

namespace smetanin_d_hoare_even_odd_batchelor {

class SmetaninDHoarSortOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit SmetaninDHoarSortOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace smetanin_d_hoare_even_odd_batchelor
