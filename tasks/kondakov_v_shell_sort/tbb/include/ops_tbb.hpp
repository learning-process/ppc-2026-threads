#pragma once

#include "kondakov_v_shell_sort/common/include/common.hpp"
#include "task/include/task.hpp"

namespace kondakov_v_shell_sort {

class KondakovVShellSortTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit KondakovVShellSortTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace kondakov_v_shell_sort
