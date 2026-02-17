#pragma once

#include "dergachev_a_graham_scan/common/include/common.hpp"
#include "task/include/task.hpp"

namespace dergachev_a_graham_scan {

class DergachevAGrahamScanSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit DergachevAGrahamScanSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace dergachev_a_graham_scan
