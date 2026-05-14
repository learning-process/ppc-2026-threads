#pragma once

#include "gaivoronskiy_m_marking_binary_components/common/include/common.hpp"
#include "task/include/task.hpp"

namespace gaivoronskiy_m_marking_binary_components {

class GaivoronskiyMMarkingBinaryComponentsSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit GaivoronskiyMMarkingBinaryComponentsSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace gaivoronskiy_m_marking_binary_components
