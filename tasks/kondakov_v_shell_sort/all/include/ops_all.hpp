#pragma once

#include "kondakov_v_shell_sort/common/include/common.hpp"
#include "task/include/task.hpp"

namespace kondakov_v_shell_sort {

class KondakovVShellSortALL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kALL;
  }
  explicit KondakovVShellSortALL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace kondakov_v_shell_sort
