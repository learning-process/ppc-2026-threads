#pragma once

#include "task/include/task.hpp"
#include "yurkin_g_graham_scan/common/include/common.hpp"

namespace yurkin_g_graham_scan {

class YurkinGGrahamScanOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit YurkinGGrahamScanOMP(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace yurkin_g_graham_scan
