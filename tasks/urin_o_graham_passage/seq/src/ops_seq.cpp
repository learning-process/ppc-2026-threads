#include "urin_o_graham_passage/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace urin_o_graham_passage {

UrinOGrahamPassageSEQ::UrinOGrahamPassageSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  // Используем не-const версию GetInput() в конструкторе
  GetInput() = in;
  GetOutput() = OutType();
}

bool UrinOGrahamPassageSEQ::ValidationImpl() {
  const auto &points = GetInput();  // В ValidationImpl мы можем использовать не-const версию
  return points.size() >= 3;
}

bool UrinOGrahamPassageSEQ::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}

// Вспомогательные статические функции
Point UrinOGrahamPassageSEQ::FindLowestPoint(const InType &points) {
  Point lowest = points[0];

  for (size_t i = 1; i < points.size(); i++) {
    if (points[i].y < lowest.y || (std::abs(points[i].y - lowest.y) < 1e-10 && points[i].x < lowest.x)) {
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
  double val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);

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
  const InType &points = GetInput();  // Получаем ссылку на входные данные

  if (points.size() < 3) {
    return false;
  }

  // Шаг 1: Находим самую нижнюю левую точку
  Point p0 = FindLowestPoint(points);

  // Шаг 2: Копируем все точки, кроме p0, для сортировки
  std::vector<Point> sorted_points;
  for (const auto &point : points) {
    if (!(point == p0)) {
      sorted_points.push_back(point);
    }
  }

  // Сортируем по полярному углу
  std::sort(sorted_points.begin(), sorted_points.end(), [this, &p0](const Point &a, const Point &b) {
    double angle_a = PolarAngle(p0, a);
    double angle_b = PolarAngle(p0, b);

    if (std::abs(angle_a - angle_b) < 1e-10) {
      return DistanceSquared(p0, a) < DistanceSquared(p0, b);
    }
    return angle_a < angle_b;
  });

  // Удаляем коллинеарные точки
  std::vector<Point> unique_angles;
  for (size_t i = 0; i < sorted_points.size(); i++) {
    if (i > 0 && std::abs(PolarAngle(p0, sorted_points[i]) - PolarAngle(p0, sorted_points[i - 1])) < 1e-10) {
      continue;
    }
    unique_angles.push_back(sorted_points[i]);
  }

  if (unique_angles.size() < 2) {
    return false;
  }

  unique_angles.insert(unique_angles.begin(), p0);

  // Шаг 3: Строим выпуклую оболочку
  std::vector<Point> hull;
  hull.push_back(unique_angles[0]);
  hull.push_back(unique_angles[1]);

  for (size_t i = 2; i < unique_angles.size(); i++) {
    while (hull.size() >= 2 && Orientation(hull[hull.size() - 2], hull.back(), unique_angles[i]) != 1) {
      hull.pop_back();
    }
    hull.push_back(unique_angles[i]);
  }

  // Сохраняем результат через не-const версию GetOutput()
  GetOutput() = hull;

  return hull.size() >= 3;
}

bool UrinOGrahamPassageSEQ::PostProcessingImpl() {
  return !GetOutput().empty();
}

}  // namespace urin_o_graham_passage
