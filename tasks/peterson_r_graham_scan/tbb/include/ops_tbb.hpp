#pragma once

#include <cstddef>

#include "peterson_r_graham_scan/common/include/common.hpp"
#include "task/include/task.hpp"

namespace peterson_r_graham_scan {

class PetersonGrahamScannerTbb : public TaskBase {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kTBB;
  }

  explicit PetersonGrahamScannerTbb(const InputValue &in);

  void LoadPoints(const PointSet &points);
  [[nodiscard]] PointSet GetConvexHull() const;

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  PointSet input_points_;
  PointSet hull_points_;
  bool external_data_provided_ = false;

  static double ComputeOrientation(const Point2D &origin, const Point2D &a, const Point2D &b);
  static double ComputeDistanceSq(const Point2D &p1, const Point2D &p2);
  static bool AreAllPointsIdentical(const PointSet &points);
  static std::size_t FindLowestPointParallel(const PointSet &points);
  static void SortPointsByAngleParallel(PointSet &points);
};

}  // namespace peterson_r_graham_scan
