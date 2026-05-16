#pragma once

#include <cstddef>
#include <utility>
#include <vector>

#include "perepelkin_i_convex_hull_graham_scan/common/include/common.hpp"
#include "task/include/task.hpp"

namespace perepelkin_i_convex_hull_graham_scan {

class PerepelkinIConvexHullGrahamScanTBB : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }
  explicit PerepelkinIConvexHullGrahamScanTBB(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static size_t FindPivotParallel(const std::vector<std::pair<double, double>> &pts);
  static void ParallelSort(std::vector<std::pair<double, double>> &data, const std::pair<double, double> &pivot);

  static void HullConstruction(std::vector<std::pair<double, double>> &hull,
                               const std::vector<std::pair<double, double>> &pts,
                               const std::pair<double, double> &pivot);
  static bool AngleCmp(const std::pair<double, double> &a, const std::pair<double, double> &b,
                       const std::pair<double, double> &pivot);
  static double Orientation(const std::pair<double, double> &p, const std::pair<double, double> &q,
                            const std::pair<double, double> &r);
};

}  // namespace perepelkin_i_convex_hull_graham_scan
