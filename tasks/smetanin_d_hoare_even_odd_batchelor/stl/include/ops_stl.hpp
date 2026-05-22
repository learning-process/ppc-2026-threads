#pragma once

#include "smetanin_d_hoare_even_odd_batchelor/common/include/common.hpp"
#include "task/include/task.hpp"

namespace smetanin_d_hoare_even_odd_batchelor {

class SmetaninDHoarSortSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit SmetaninDHoarSortSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace smetanin_d_hoare_even_odd_batchelor
