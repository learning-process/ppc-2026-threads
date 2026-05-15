#include "kosolapov_v_calc_mult_integrals_m_rectangles/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <execution>
#include <numeric>
#include <ranges>
#include <tuple>

#include "kosolapov_v_calc_mult_integrals_m_rectangles/common/include/common.hpp"

namespace kosolapov_v_calc_mult_integrals_m_rectangles {

KosolapovVCalcMultIntegralsMRectanglesSTL::KosolapovVCalcMultIntegralsMRectanglesSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = InType(in);
  GetOutput() = 0.0;
}

bool KosolapovVCalcMultIntegralsMRectanglesSTL::ValidationImpl() {
  int steps = std::get<0>(GetInput());
  int func_id = std::get<1>(GetInput());
  return steps > 0 && func_id >= 1 && func_id <= 4;
}

bool KosolapovVCalcMultIntegralsMRectanglesSTL::PreProcessingImpl() {
  return true;
}

bool KosolapovVCalcMultIntegralsMRectanglesSTL::RunImpl() {
  int steps = std::get<0>(GetInput());
  int func_id = std::get<1>(GetInput());
  std::tuple<double, double, double, double> temp = GetBounds(func_id);
  double a = std::get<0>(temp);
  double b = std::get<1>(temp);
  double c = std::get<2>(temp);
  double d = std::get<3>(temp);
  double integral = RectanglesIntegral(func_id, steps, a, b, c, d);
  GetOutput() = integral;
  return true;
}

bool KosolapovVCalcMultIntegralsMRectanglesSTL::PostProcessingImpl() {
  return true;
}

double KosolapovVCalcMultIntegralsMRectanglesSTL::Function1(double x, double y) {
  // f(x,y) = x^2 + y^2
  return (x * x) + (y * y);
}
double KosolapovVCalcMultIntegralsMRectanglesSTL::Function2(double x, double y) {
  // f(x,y) = sin(x) * cos(y)
  return std::sin(x) * std::cos(y);
}
double KosolapovVCalcMultIntegralsMRectanglesSTL::Function3(double x, double y) {
  // f(x,y) = exp(-(x^2 + y^2))
  return std::exp(-((x * x) + (y * y)));
}
double KosolapovVCalcMultIntegralsMRectanglesSTL::Function4(double x, double y) {
  // f(x,y) = sin(x + y)
  return std::sin(x + y);
}
double KosolapovVCalcMultIntegralsMRectanglesSTL::CallFunction(int func_id, double x, double y) {
  switch (func_id) {
    case 1:
      return Function1(x, y);
    case 2:
      return Function2(x, y);
    case 3:
      return Function3(x, y);
    case 4:
      return Function4(x, y);
    default:
      return Function1(x, y);
  }
}
std::tuple<double, double, double, double> KosolapovVCalcMultIntegralsMRectanglesSTL::GetBounds(int func_id) {
  switch (func_id) {
    case 1:
      return {0.0, 1.0, 0.0, 1.0};
    case 2:
      return {0.0, kPi, 0.0, kPi / 2.0};
    case 3:
      return {-1.0, 1.0, -1.0, 1.0};
    case 4:
      return {0.0, kPi, 0.0, kPi};
    default:
      return {0.0, 1.0, 0.0, 1.0};
  }
}
double KosolapovVCalcMultIntegralsMRectanglesSTL::RectanglesIntegral(int func_id, int steps, double a, double b,
                                                                     double c, double d) {
  double hx = (b - a) / steps;
  double hy = (d - c) / steps;
  size_t total = static_cast<size_t>(steps) * steps;

  auto indices = std::views::iota(0ULL, total);
  double sum =
      std::transform_reduce(std::execution::par, indices.begin(), indices.end(), 0.0, std::plus<>(), [&](size_t idx) {
    int i = static_cast<int>(idx / steps);
    int j = static_cast<int>(idx % steps);
    double x = a + ((i + 0.5) * hx);
    double y = c + ((j + 0.5) * hy);
    return CallFunction(func_id, x, y);
  });

  return sum * hx * hy;
}

}  // namespace kosolapov_v_calc_mult_integrals_m_rectangles
