#pragma once

#include <utility>
#include <vector>

#include "dergachev_a_graham_scan/common/include/common.hpp"
#include "task/include/task.hpp"

namespace dergachev_a_graham_scan {

class DergachevAGrahamScanSTL : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSTL;
  }
  explicit DergachevAGrahamScanSTL(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  std::vector<std::pair<double, double>> points_;
  std::vector<std::pair<double, double>> hull_;
};

}  // namespace dergachev_a_graham_scan
