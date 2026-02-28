#pragma once

#include <vector>

#include "peterson_r_graham_scan_omp/common/include/common.hpp"
#include "task/include/task.hpp"

namespace peterson_r_graham_scan_omp {

struct Point {
  double x;
  double y;
};

class PetersonRGrahamScanOMP : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kOMP;
  }
  explicit PetersonRGrahamScanOMP(const InType &in);

  void SetPoints(const std::vector<Point> &pts);
  [[nodiscard]] std::vector<Point> GetHull() const;

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  std::vector<Point> input_points_;
  std::vector<Point> convex_hull_;
  bool is_custom_input_ = false;
};

}  // namespace peterson_r_graham_scan_omp
