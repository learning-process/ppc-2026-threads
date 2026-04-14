#pragma once

#include <cstdint>
#include <vector>

#include "rychkova_gauss/common/include/common.hpp"
#include "task/include/task.hpp"

namespace rychkova_gauss {

class RychkovaGaussALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit RychkovaGaussALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
  std::vector<std::uint8_t> loutput_;  // результаты лок. вычислений
};

}  // namespace rychkova_gauss
