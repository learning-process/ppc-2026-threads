#pragma once

#include <cstddef>

#include "kutuzov_i_convex_hull_jarvis/common/include/common.hpp"
#include "task/include/task.hpp"

namespace kutuzov_i_convex_hull_jarvis {

class KutuzovITestConvexHullSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit KutuzovITestConvexHullSEQ(const InType &in);

 private:
  static double DistanceSquared(double a_x, double a_y, double b_x, double b_y);
  static double CrossProduct(double o_x, double o_y, double a_x, double a_y, double b_x, double b_y);

  static size_t FindLeftmostPoint(const InType &input);
  static bool IsBetterPoint(double cross, double epsilon, double current_x, double current_y, double i_x, double i_y,
                            double next_x, double next_y);

  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace kutuzov_i_convex_hull_jarvis
