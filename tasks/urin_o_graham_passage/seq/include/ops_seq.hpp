#pragma once

#include <vector>

#include "task/include/task.hpp"
#include "urin_o_graham_passage/common/include/common.hpp"

namespace urin_o_graham_passage {

class UrinOGrahamPassageSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }

  explicit UrinOGrahamPassageSEQ(const InType &in);

  // Статические вспомогательные функции
  [[nodiscard]] static Point FindLowestPoint(const InType &points);
  [[nodiscard]] static double PolarAngle(const Point &base, const Point &p);
  [[nodiscard]] static int Orientation(const Point &p, const Point &q, const Point &r);
  [[nodiscard]] static double DistanceSquared(const Point &p1, const Point &p2);

 protected:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace urin_o_graham_passage
