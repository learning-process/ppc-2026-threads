#pragma once

#include "konstantinov_s_graham/common/include/common.hpp"
#include "task/include/task.hpp"

namespace konstantinov_a_graham {

class KonstantinovAGrahamSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit KonstantinovAGrahamSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace konstantinov_a_graham
