#pragma once

#include "task/include/task.hpp"

namespace ppc::util::detail {

constexpr bool IsCompiledTaskImplEnabled(ppc::task::TypeOfTask task_type) {
#ifndef PPC_TASK_IMPL_FILTERED
  (void)task_type;
  return true;
#else
  switch (task_type) {
    case ppc::task::TypeOfTask::kALL:
#  ifdef PPC_TASK_IMPL_ALL
      return true;
#  else
      return false;
#  endif
    case ppc::task::TypeOfTask::kMPI:
#  ifdef PPC_TASK_IMPL_MPI
      return true;
#  else
      return false;
#  endif
    case ppc::task::TypeOfTask::kOMP:
#  ifdef PPC_TASK_IMPL_OMP
      return true;
#  else
      return false;
#  endif
    case ppc::task::TypeOfTask::kSEQ:
#  ifdef PPC_TASK_IMPL_SEQ
      return true;
#  else
      return false;
#  endif
    case ppc::task::TypeOfTask::kSTL:
#  ifdef PPC_TASK_IMPL_STL
      return true;
#  else
      return false;
#  endif
    case ppc::task::TypeOfTask::kTBB:
#  ifdef PPC_TASK_IMPL_TBB
      return true;
#  else
      return false;
#  endif
    case ppc::task::TypeOfTask::kUnknown:
      return true;
  }
  return true;
#endif
}

template <ppc::task::TypeOfTask TaskType>
constexpr bool IsCompiledTaskImplEnabled() {
  return IsCompiledTaskImplEnabled(TaskType);
}

}  // namespace ppc::util::detail
