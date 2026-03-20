#pragma once

#include <utility>

#include "perepelkin_i_convex_hull_graham_scan/common/include/common.hpp"
#include "task/include/task.hpp"

namespace perepelkin_i_convex_hull_graham_scan {

class PerepelkinIConvexHullGrahamScanSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit PerepelkinIConvexHullGrahamScanSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static bool AngleCmp(const std::pair<double, double> &a, const std::pair<double, double> &b,
                       const std::pair<double, double> &pivot);
  static double Orientation(const std::pair<double, double> &p, const std::pair<double, double> &q,
                            const std::pair<double, double> &r);
};

}  // namespace perepelkin_i_convex_hull_graham_scan
