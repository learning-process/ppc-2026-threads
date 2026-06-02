#pragma once

#include "task/include/task.hpp"

namespace ppc::util::detail {

#ifdef PPC_TASK_IMPL_FILTERED

#  ifdef PPC_TASK_IMPL_ALL
constexpr bool kTaskImplAllEnabled = true;
#  else
constexpr bool kTaskImplAllEnabled = false;
#  endif

#  ifdef PPC_TASK_IMPL_MPI
constexpr bool kTaskImplMpiEnabled = true;
#  else
constexpr bool kTaskImplMpiEnabled = false;
#  endif

#  ifdef PPC_TASK_IMPL_OMP
constexpr bool kTaskImplOmpEnabled = true;
#  else
constexpr bool kTaskImplOmpEnabled = false;
#  endif

#  ifdef PPC_TASK_IMPL_SEQ
constexpr bool kTaskImplSeqEnabled = true;
#  else
constexpr bool kTaskImplSeqEnabled = false;
#  endif

#  ifdef PPC_TASK_IMPL_STL
constexpr bool kTaskImplStlEnabled = true;
#  else
constexpr bool kTaskImplStlEnabled = false;
#  endif

#  ifdef PPC_TASK_IMPL_TBB
constexpr bool kTaskImplTbbEnabled = true;
#  else
constexpr bool kTaskImplTbbEnabled = false;
#  endif

#endif

constexpr bool IsCompiledTaskImplEnabled(ppc::task::TypeOfTask task_type) {
#ifndef PPC_TASK_IMPL_FILTERED
  (void)task_type;
  return true;
#else
  switch (task_type) {
    case ppc::task::TypeOfTask::kALL:
      return kTaskImplAllEnabled;
    case ppc::task::TypeOfTask::kMPI:
      return kTaskImplMpiEnabled;
    case ppc::task::TypeOfTask::kOMP:
      return kTaskImplOmpEnabled;
    case ppc::task::TypeOfTask::kSEQ:
      return kTaskImplSeqEnabled;
    case ppc::task::TypeOfTask::kSTL:
      return kTaskImplStlEnabled;
    case ppc::task::TypeOfTask::kTBB:
      return kTaskImplTbbEnabled;
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
