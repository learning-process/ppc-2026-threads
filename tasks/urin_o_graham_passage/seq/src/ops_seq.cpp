#include "urin_o_graham_passage/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace urin_o_graham_passage {

UrinOGrahamPassageSEQ::UrinOGrahamPassageSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType();
}

bool UrinOGrahamPassageSEQ::ValidationImpl() {
  const auto &points = GetInput();

  // Проверяем, что точек достаточно
  if (points.size() < 3) {
    return false;
  }

  // Проверяем, что не все точки одинаковые
  const Point &first = points[0];
  bool all_same = true;
  for (size_t i = 1; i < points.size(); i++) {
    if (!(points[i] == first)) {
      all_same = false;
      break;
    }
  }

  // Если все точки одинаковые - валидация не проходит
  if (all_same) {
    return false;
  }

  return true;
}

bool UrinOGrahamPassageSEQ::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}

Point UrinOGrahamPassageSEQ::FindLowestPoint(const InType &points) {
  Point lowest = points[0];

  for (size_t i = 1; i < points.size(); i++) {
    if (points[i].y < lowest.y - 1e-10 ||
        (std::abs(points[i].y - lowest.y) < 1e-10 && points[i].x < lowest.x - 1e-10)) {
      lowest = points[i];
    }
  }

  return lowest;
}

double UrinOGrahamPassageSEQ::PolarAngle(const Point &base, const Point &p) {
  double dx = p.x - base.x;
  double dy = p.y - base.y;

  if (std::abs(dx) < 1e-10 && std::abs(dy) < 1e-10) {
    return -1e10;
  }

  return std::atan2(dy, dx);
}

int UrinOGrahamPassageSEQ::Orientation(const Point &p, const Point &q, const Point &r) {
  // Правильная формула для определения поворота
  double val = (q.x - p.x) * (r.y - p.y) - (q.y - p.y) * (r.x - p.x);

  if (std::abs(val) < 1e-10) {
    return 0;
  }
  return (val > 0) ? 1 : -1;
}

double UrinOGrahamPassageSEQ::DistanceSquared(const Point &p1, const Point &p2) {
  double dx = p2.x - p1.x;
  double dy = p2.y - p1.y;
  return dx * dx + dy * dy;
}

bool UrinOGrahamPassageSEQ::RunImpl() {
  const InType &points = GetInput();

  if (points.size() < 3) {
    return false;
  }

  // Шаг 1: Находим самую нижнюю левую точку
  Point p0 = FindLowestPoint(points);

  // Шаг 2: Копируем все точки, кроме p0
  std::vector<Point> other_points;
  for (const auto &point : points) {
    if (!(point == p0)) {
      other_points.push_back(point);
    }
  }

  // Шаг 3: Сортируем по полярному углу относительно p0
  std::sort(other_points.begin(), other_points.end(), [this, &p0](const Point &a, const Point &b) {
    double angle_a = PolarAngle(p0, a);
    double angle_b = PolarAngle(p0, b);

    if (std::abs(angle_a - angle_b) < 1e-10) {
      return DistanceSquared(p0, a) < DistanceSquared(p0, b);
    }
    return angle_a < angle_b;
  });

  // Шаг 4: Проверяем, не все ли точки коллинеарны
  bool all_collinear = true;
  for (size_t i = 1; i < other_points.size(); i++) {
    if (Orientation(p0, other_points[0], other_points[i]) != 0) {
      all_collinear = false;
      break;
    }
  }

  if (all_collinear) {
    // Все точки на одной прямой - возвращаем две крайние
    GetOutput() = {p0, other_points.back()};
    return true;
  }

  // Шаг 5: Классический алгоритм Грэхема
  std::vector<Point> hull;
  hull.push_back(p0);
  hull.push_back(other_points[0]);

  for (size_t i = 1; i < other_points.size(); i++) {
    while (hull.size() >= 2) {
      Point &p = hull[hull.size() - 2];
      Point &q = hull[hull.size() - 1];
      int orient = Orientation(p, q, other_points[i]);

      if (orient <= 0) {  // Правый поворот или коллинеарны - удаляем
        hull.pop_back();
      } else {
        break;
      }
    }
    hull.push_back(other_points[i]);
  }

  GetOutput() = hull;
  return true;
}

bool UrinOGrahamPassageSEQ::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace urin_o_graham_passage
