#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "urin_o_graham_passage/common/include/common.hpp"
#include "urin_o_graham_passage/seq/include/ops_seq.hpp"

namespace urin_o_graham_passage {

// Вспомогательная функция для проверки выпуклости оболочки
bool IsConvexHull(const std::vector<Point> &hull) {
  if (hull.size() < 3) {
    return true;  // Для 1-2 точек считаем корректным
  }

  for (size_t i = 0; i < hull.size(); i++) {
    size_t prev = (i == 0) ? hull.size() - 1 : i - 1;
    size_t next = (i + 1) % hull.size();

    if (UrinOGrahamPassageSEQ::Orientation(hull[prev], hull[i], hull[next]) < 0) {
      return false;
    }
  }
  return true;
}

// Тест с пустым входом
TEST(UrinOGrahamPassageSeq, EmptyInput) {
  InType empty_points;
  auto task = std::make_shared<UrinOGrahamPassageSEQ>(empty_points);

  EXPECT_FALSE(task->Validation());
  EXPECT_TRUE(task->GetOutput().empty());
}

// Тест с одной точкой
TEST(UrinOGrahamPassageSeq, SinglePoint) {
  InType pts = {Point(5.0, 3.0)};
  auto task = std::make_shared<UrinOGrahamPassageSEQ>(pts);

  EXPECT_FALSE(task->Validation());
  EXPECT_TRUE(task->GetOutput().empty());
}

// Тест с двумя точками
TEST(UrinOGrahamPassageSeq, TwoDistinctPoints) {
  InType pts = {Point(0.0, 0.0), Point(3.0, 4.0)};
  auto task = std::make_shared<UrinOGrahamPassageSEQ>(pts);

  EXPECT_FALSE(task->Validation());
  EXPECT_TRUE(task->GetOutput().empty());
}

// Тест с коллинеарными точками на прямой
TEST(UrinOGrahamPassageSeq, CollinearPoints) {
  InType pts = {Point(0.0, 0.0), Point(1.0, 0.0), Point(2.0, 0.0), Point(3.0, 0.0), Point(4.0, 0.0)};
  auto task = std::make_shared<UrinOGrahamPassageSEQ>(pts);

  EXPECT_TRUE(task->Validation());
  EXPECT_TRUE(task->PreProcessing());
  EXPECT_TRUE(task->Run());
  EXPECT_TRUE(task->PostProcessing());

  const auto &hull = task->GetOutput();
  // ИСПРАВЛЕНО: сравниваем size_t с size_t
  EXPECT_EQ(hull.size(), size_t(2));  // Должны получить две крайние точки
  EXPECT_TRUE(IsConvexHull(hull));
}

// Тест с треугольником
TEST(UrinOGrahamPassageSeq, TrianglePoints) {
  InType pts = {Point(0.0, 0.0), Point(4.0, 0.0), Point(2.0, 3.0)};
  auto task = std::make_shared<UrinOGrahamPassageSEQ>(pts);

  EXPECT_TRUE(task->Validation());
  EXPECT_TRUE(task->PreProcessing());
  EXPECT_TRUE(task->Run());
  EXPECT_TRUE(task->PostProcessing());

  const auto &hull = task->GetOutput();
  EXPECT_EQ(hull.size(), size_t(3));
  EXPECT_TRUE(IsConvexHull(hull));
}

// Тест с квадратом
TEST(UrinOGrahamPassageSeq, SquarePoints) {
  InType pts = {Point(0.0, 0.0), Point(4.0, 0.0), Point(4.0, 4.0), Point(0.0, 4.0)};
  auto task = std::make_shared<UrinOGrahamPassageSEQ>(pts);

  EXPECT_TRUE(task->Validation());
  EXPECT_TRUE(task->PreProcessing());
  EXPECT_TRUE(task->Run());
  EXPECT_TRUE(task->PostProcessing());

  const auto &hull = task->GetOutput();
  EXPECT_EQ(hull.size(), size_t(4));
  EXPECT_TRUE(IsConvexHull(hull));
}

// Тест с квадратом и внутренней точкой
TEST(UrinOGrahamPassageSeq, SquareWithInteriorPoint) {
  InType pts = {Point(0.0, 0.0), Point(4.0, 0.0), Point(4.0, 4.0), Point(0.0, 4.0), Point(2.0, 2.0)};
  auto task = std::make_shared<UrinOGrahamPassageSEQ>(pts);

  EXPECT_TRUE(task->Validation());
  EXPECT_TRUE(task->PreProcessing());
  EXPECT_TRUE(task->Run());
  EXPECT_TRUE(task->PostProcessing());

  const auto &hull = task->GetOutput();
  EXPECT_EQ(hull.size(), size_t(4));
  EXPECT_TRUE(IsConvexHull(hull));
}

// Тест с прямоугольником и коллинеарными точками на сторонах
TEST(UrinOGrahamPassageSeq, RectangleWithCollinearPoints) {
  InType pts = {Point(0.0, 0.0), Point(1.0, 0.0), Point(2.0, 0.0), Point(3.0, 0.0),
                Point(3.0, 1.0), Point(2.0, 1.0), Point(1.0, 1.0), Point(0.0, 1.0)};
  auto task = std::make_shared<UrinOGrahamPassageSEQ>(pts);

  EXPECT_TRUE(task->Validation());
  EXPECT_TRUE(task->PreProcessing());
  EXPECT_TRUE(task->Run());
  EXPECT_TRUE(task->PostProcessing());

  const auto &hull = task->GetOutput();
  EXPECT_EQ(hull.size(), size_t(4));  // Должны получить 4 вершины прямоугольника
  EXPECT_TRUE(IsConvexHull(hull));
}

// Тест с повторяющимися точками
TEST(UrinOGrahamPassageSeq, AllIdenticalPoints) {
  InType pts = {Point(3.0, 3.0), Point(3.0, 3.0), Point(3.0, 3.0), Point(3.0, 3.0), Point(3.0, 3.0)};
  auto task = std::make_shared<UrinOGrahamPassageSEQ>(pts);

  EXPECT_FALSE(task->Validation());  // Валидация должна вернуть false
  EXPECT_TRUE(task->GetOutput().empty());
}

// Тест с точками на границе
TEST(UrinOGrahamPassageSeq, PointOnBoundary) {
  InType pts = {Point(0.0, 0.0), Point(4.0, 0.0), Point(2.0, 0.0), Point(4.0, 4.0), Point(0.0, 4.0)};
  auto task = std::make_shared<UrinOGrahamPassageSEQ>(pts);

  EXPECT_TRUE(task->Validation());
  EXPECT_TRUE(task->PreProcessing());
  EXPECT_TRUE(task->Run());
  EXPECT_TRUE(task->PostProcessing());

  const auto &hull = task->GetOutput();
  EXPECT_EQ(hull.size(), size_t(4));  // Точка (2,0) лежит на стороне
  EXPECT_TRUE(IsConvexHull(hull));
}

// Тест с вертикальной линией
TEST(UrinOGrahamPassageSeq, VerticalCollinear) {
  InType pts = {Point(0.0, 0.0), Point(0.0, 1.0), Point(0.0, 2.0), Point(0.0, 5.0)};
  auto task = std::make_shared<UrinOGrahamPassageSEQ>(pts);

  EXPECT_TRUE(task->Validation());
  EXPECT_TRUE(task->PreProcessing());
  EXPECT_TRUE(task->Run());
  EXPECT_TRUE(task->PostProcessing());

  const auto &hull = task->GetOutput();
  EXPECT_EQ(hull.size(), size_t(2));  // Две крайние точки
  EXPECT_TRUE(IsConvexHull(hull));
}

// Тест с большим количеством случайных точек
TEST(UrinOGrahamPassageSeq, LargeRandomSet) {
  InType pts;
  const int num_points = 100;

  // Генерируем точки на окружности
  for (int i = 0; i < num_points; ++i) {
    double angle = 2.0 * 3.14159 * i / num_points;
    pts.emplace_back(cos(angle) * 10.0, sin(angle) * 10.0);
  }

  auto task = std::make_shared<UrinOGrahamPassageSEQ>(pts);

  EXPECT_TRUE(task->Validation());
  EXPECT_TRUE(task->PreProcessing());
  EXPECT_TRUE(task->Run());
  EXPECT_TRUE(task->PostProcessing());

  const auto &hull = task->GetOutput();
  // ИСПРАВЛЕНО: сравниваем беззнаковые типы
  EXPECT_GE(hull.size(), size_t(3));
  EXPECT_TRUE(IsConvexHull(hull));
}

// Тест с шестиугольником и центром
TEST(UrinOGrahamPassageSeq, HexagonWithCenter) {
  InType pts = {Point(2.0, 0.0),    Point(1.0, 1.73),  Point(-1.0, 1.73), Point(-2.0, 0.0),
                Point(-1.0, -1.73), Point(1.0, -1.73), Point(0.0, 0.0)};
  auto task = std::make_shared<UrinOGrahamPassageSEQ>(pts);

  EXPECT_TRUE(task->Validation());
  EXPECT_TRUE(task->PreProcessing());
  EXPECT_TRUE(task->Run());
  EXPECT_TRUE(task->PostProcessing());

  const auto &hull = task->GetOutput();
  EXPECT_EQ(hull.size(), size_t(6));  // Шестиугольник без центра
  EXPECT_TRUE(IsConvexHull(hull));
}

}  // namespace urin_o_graham_passage
