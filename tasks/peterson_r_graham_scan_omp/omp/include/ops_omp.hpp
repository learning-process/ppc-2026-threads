#pragma once

#include <cstddef>
#include <vector>

#include "peterson_r_graham_scan_omp/common/include/common.hpp"
#include "task/include/task.hpp"

namespace peterson_r_graham_scan_omp {

struct Point2D {
  double coord_x;
  double coord_y;

  Point2D() = default;
  Point2D(double x_val, double y_val) : coord_x(x_val), coord_y(y_val) {}
};

using PointSet = std::vector<Point2D>;

class PetersonGrahamScannerOMP : public ppc::task::Task<InputValue, OutputValue> {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }

  explicit PetersonGrahamScannerOMP(const InputValue &in);

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
  static void ParallelMergeSort(PointSet &points, int left, int right, const Point2D &origin);
};

}  // namespace peterson_r_graham_scan_omp
