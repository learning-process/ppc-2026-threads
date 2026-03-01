#pragma once

#include <cstddef>
#include <vector>

#include "peterson_r_graham_scan_seq/common/include/common.hpp"
#include "task/include/task.hpp"

namespace peterson_r_graham_scan_seq {

struct Coordinate2D {
  double x;
  double y;

  Coordinate2D() = default;
  Coordinate2D(double x_val, double y_val) : x(x_val), y(y_val) {}
};

using PointCloud = std::vector<Coordinate2D>;

class PetersonRGrahamScanSeq : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSeq;
  }

  explicit PetersonRGrahamScanSeq(const InType &in);

  void SetPoints(const PointCloud &cloud);
  [[nodiscard]] PointCloud GetHull() const;

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  PointCloud dataset_;
  PointCloud boundary_;
  bool user_provided_data_ = false;

  static double CrossProduct(const Coordinate2D &origin, const Coordinate2D &a, const Coordinate2D &b);
  static double SquaredDistance(const Coordinate2D &p1, const Coordinate2D &p2);
  static bool IsUniformCloud(const PointCloud &cloud);
  static std::size_t FindAnchorIndex(const PointCloud &cloud);
  static void ArrangeByPolarAngle(PointCloud &cloud);
};

}  // namespace peterson_r_graham_scan_seq
