#pragma once

#include <tuple>

#include "kosolapov_v_calc_mult_integrals_m_rectangles/common/include/common.hpp"
#include "task/include/task.hpp"

namespace kosolapov_v_calc_mult_integrals_m_rectangles {

class KosolapovVCalcMultIntegralsMRectanglesSEQ : public BaseTask {
 public:
  static constexpr ppc::task::TypeOfTask GetStaticTypeOfTask() {
    return ppc::task::TypeOfTask::kSEQ;
  }
  explicit KosolapovVCalcMultIntegralsMRectanglesSEQ(const InType &in);

 private:
  bool ValidationImpl() override;
  bool PreProcessingImpl() override;
  bool RunImpl() override;
  bool PostProcessingImpl() override;

  static double Function1(double x, double y);
  static double Function2(double x, double y);
  static double Function3(double x, double y);
  static double Function4(double x, double y);
  static double CallFunction(int func_id, double x, double y);
  static std::tuple<double, double, double, double> GetBounds(int func_id);
  static double RectanglesIntegral(int func_id, int steps, double a, double b, double c, double d);
};

}  // namespace kosolapov_v_calc_mult_integrals_m_rectangles
