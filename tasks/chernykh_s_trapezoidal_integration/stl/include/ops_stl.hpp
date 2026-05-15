#pragma once
#include <cstddef>
#include <vector>

#include "chernykh_s_trapezoidal_integration/common/include/common.hpp"
#include "task/include/task.hpp"
namespace chernykh_s_trapezoidal_integration {

class ChernykhSTrapezoidalIntegrationSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit ChernykhSTrapezoidalIntegrationSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace chernykh_s_trapezoidal_integration
