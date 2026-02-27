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

  // Публичные вспомогательные функции для тестов
  static Point FindLowestPoint(const InType &points);
  static double PolarAngle(const Point &base, const Point &p);
  static int Orientation(const Point &p, const Point &q, const Point &r);
  static double DistanceSquared(const Point &p1, const Point &p2);

  // Константные геттеры для тестов
  const InType &GetInputPoints() const {
    return const_cast<UrinOGrahamPassageSEQ *>(this)->GetInput();
  }

  const OutType &GetHull() const {
    return const_cast<UrinOGrahamPassageSEQ *>(this)->GetOutput();
  }

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;
};

}  // namespace urin_o_graham_passage
