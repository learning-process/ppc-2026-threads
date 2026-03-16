#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <memory>
#include <numbers>
#include <vector>

#include "urin_o_graham_passage/common/include/common.hpp"
#include "urin_o_graham_passage/omp/include/ops_omp.hpp"  // ИЗМЕНЕНО: вместо seq

namespace urin_o_graham_passage {
namespace {

bool IsConvexHull(const std::vector<Point> &hull) {
  if (hull.size() < 3) {
    return true;
  }

  for (size_t i = 0; i < hull.size(); ++i) {
    size_t prev = (i == 0) ? hull.size() - 1 : i - 1;
    size_t next = (i + 1) % hull.size();

    if (UrinOGrahamPassageOMP::Orientation(hull[prev], hull[i], hull[next]) < 0) {  // ИЗМЕНЕНО
      return false;
    }
  }
  return true;
}

// Еще более мелкие функции для каждого шага
bool ValidateTask(const std::shared_ptr<UrinOGrahamPassageOMP> &task) {  // ИЗМЕНЕНО
  return task->Validation();
}

bool PreProcessTask(const std::shared_ptr<UrinOGrahamPassageOMP> &task) {  // ИЗМЕНЕНО
  return task->PreProcessing();
}

bool RunTask(const std::shared_ptr<UrinOGrahamPassageOMP> &task) {  // ИЗМЕНЕНО
  return task->Run();
}

bool PostProcessTask(const std::shared_ptr<UrinOGrahamPassageOMP> &task) {  // ИЗМЕНЕНО
  return task->PostProcessing();
}

// Отдельные функции для проверок с EXPECT_
void ExpectValidation(const std::shared_ptr<UrinOGrahamPassageOMP> &task) {  // ИЗМЕНЕНО
  EXPECT_TRUE(ValidateTask(task));
}

void ExpectPreProcessing(const std::shared_ptr<UrinOGrahamPassageOMP> &task) {  // ИЗМЕНЕНО
  EXPECT_TRUE(PreProcessTask(task));
}

void ExpectRun(const std::shared_ptr<UrinOGrahamPassageOMP> &task) {  // ИЗМЕНЕНО
  EXPECT_TRUE(RunTask(task));
}

void ExpectPostProcessing(const std::shared_ptr<UrinOGrahamPassageOMP> &task) {  // ИЗМЕНЕНО
  EXPECT_TRUE(PostProcessTask(task));
}

void ExecuteTaskPipeline(const std::shared_ptr<UrinOGrahamPassageOMP> &task) {  // ИЗМЕНЕНО
  ExpectValidation(task);
  ExpectPreProcessing(task);
  ExpectRun(task);
  ExpectPostProcessing(task);
}

void CheckHullSize(const std::vector<Point> &hull, size_t expected_size) {
  EXPECT_EQ(hull.size(), expected_size);
}

void CheckHullConvexity(const std::vector<Point> &hull) {
  EXPECT_TRUE(IsConvexHull(hull));
}

void VerifyHull(const std::vector<Point> &hull, size_t expected_size) {
  CheckHullSize(hull, expected_size);
  CheckHullConvexity(hull);
}

void RunAndCheckHull(const InType &points, size_t expected_size) {
  auto task = std::make_shared<UrinOGrahamPassageOMP>(points);  // ИЗМЕНЕНО
  ExecuteTaskPipeline(task);
  VerifyHull(task->GetOutput(), expected_size);
}

void RunAndExpectFailure(const InType &points) {
  auto task = std::make_shared<UrinOGrahamPassageOMP>(points);  // ИЗМЕНЕНО
  EXPECT_FALSE(task->Validation());
  EXPECT_TRUE(task->GetOutput().empty());
}

// Тесты
TEST(UrinOGrahamPassageOmp, EmptyInput) {  // ИЗМЕНЕНО: Seq -> Omp
  RunAndExpectFailure({});
}

TEST(UrinOGrahamPassageOmp, SinglePoint) {  // ИЗМЕНЕНО
  RunAndExpectFailure({Point(5.0, 3.0)});
}

TEST(UrinOGrahamPassageOmp, TwoDistinctPoints) {  // ИЗМЕНЕНО
  RunAndExpectFailure({Point(0.0, 0.0), Point(3.0, 4.0)});
}

TEST(UrinOGrahamPassageOmp, CollinearPoints) {  // ИЗМЕНЕНО
  InType pts = {Point(0.0, 0.0), Point(1.0, 0.0), Point(2.0, 0.0), Point(3.0, 0.0), Point(4.0, 0.0)};
  RunAndCheckHull(pts, 2);
}

TEST(UrinOGrahamPassageOmp, TrianglePoints) {  // ИЗМЕНЕНО
  InType pts = {Point(0.0, 0.0), Point(4.0, 0.0), Point(2.0, 3.0)};
  RunAndCheckHull(pts, 3);
}

TEST(UrinOGrahamPassageOmp, SquarePoints) {  // ИЗМЕНЕНО
  InType pts = {Point(0.0, 0.0), Point(4.0, 0.0), Point(4.0, 4.0), Point(0.0, 4.0)};
  RunAndCheckHull(pts, 4);
}

TEST(UrinOGrahamPassageOmp, SquareWithInteriorPoint) {  // ИЗМЕНЕНО
  InType pts = {Point(0.0, 0.0), Point(4.0, 0.0), Point(4.0, 4.0), Point(0.0, 4.0), Point(2.0, 2.0)};
  RunAndCheckHull(pts, 4);
}

TEST(UrinOGrahamPassageOmp, RectangleWithCollinearPoints) {  // ИЗМЕНЕНО
  InType pts = {Point(0.0, 0.0), Point(1.0, 0.0), Point(2.0, 0.0), Point(3.0, 0.0),
                Point(3.0, 1.0), Point(2.0, 1.0), Point(1.0, 1.0), Point(0.0, 1.0)};
  RunAndCheckHull(pts, 4);
}

TEST(UrinOGrahamPassageOmp, AllIdenticalPoints) {  // ИЗМЕНЕНО
  InType pts = {Point(3.0, 3.0), Point(3.0, 3.0), Point(3.0, 3.0), Point(3.0, 3.0), Point(3.0, 3.0)};
  RunAndExpectFailure(pts);
}

TEST(UrinOGrahamPassageOmp, PointOnBoundary) {  // ИЗМЕНЕНО
  InType pts = {Point(0.0, 0.0), Point(4.0, 0.0), Point(2.0, 0.0), Point(4.0, 4.0), Point(0.0, 4.0)};
  RunAndCheckHull(pts, 4);
}

TEST(UrinOGrahamPassageOmp, VerticalCollinear) {  // ИЗМЕНЕНО
  InType pts = {Point(0.0, 0.0), Point(0.0, 1.0), Point(0.0, 2.0), Point(0.0, 5.0)};
  RunAndCheckHull(pts, 2);
}

TEST(UrinOGrahamPassageOmp, LargeRandomSet) {  // ИЗМЕНЕНО
  InType pts;
  const int num_points = 100;
  pts.reserve(static_cast<size_t>(num_points));

  for (int i = 0; i < num_points; ++i) {
    const double angle = 2.0 * std::numbers::pi * static_cast<double>(i) / static_cast<double>(num_points);
    pts.emplace_back(std::cos(angle) * 10.0, std::sin(angle) * 10.0);
  }

  RunAndCheckHull(pts, 100);
}

TEST(UrinOGrahamPassageOmp, HexagonWithCenter) {  // ИЗМЕНЕНО
  InType pts = {Point(2.0, 0.0),    Point(1.0, 1.73),  Point(-1.0, 1.73), Point(-2.0, 0.0),
                Point(-1.0, -1.73), Point(1.0, -1.73), Point(0.0, 0.0)};
  RunAndCheckHull(pts, 6);
}

}  // namespace
}  // namespace urin_o_graham_passage
