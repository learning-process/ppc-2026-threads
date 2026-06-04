#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <memory>
#include <numbers>
#include <vector>

#include "urin_o_graham_passage/all/include/ops_all.hpp"
#include "urin_o_graham_passage/common/include/common.hpp"
#include "urin_o_graham_passage/omp/include/ops_omp.hpp"
#include "urin_o_graham_passage/seq/include/ops_seq.hpp"
#include "urin_o_graham_passage/stl/include/ops_stl.hpp"
#include "urin_o_graham_passage/tbb/include/ops_tbb.hpp"

namespace urin_o_graham_passage {
namespace {

double Orientation(const Point &p, const Point &q, const Point &r) {
  return ((q.x - p.x) * (r.y - p.y)) - ((q.y - p.y) * (r.x - p.x));
}

bool IsConvexHull(const std::vector<Point> &hull) {
  if (hull.size() < 3) {
    return true;
  }

  for (size_t i = 0; i < hull.size(); ++i) {
    size_t prev = (i == 0) ? hull.size() - 1 : i - 1;
    size_t next = (i + 1) % hull.size();

    if (Orientation(hull[prev], hull[i], hull[next]) < -1e-10) {
      return false;
    }
  }
  return true;
}

template <class TaskType>
bool ValidateTask(const std::shared_ptr<TaskType> &task) {
  return task->Validation();
}

template <class TaskType>
bool PreProcessTask(const std::shared_ptr<TaskType> &task) {
  return task->PreProcessing();
}

template <class TaskType>
bool RunTask(const std::shared_ptr<TaskType> &task) {
  return task->Run();
}

template <class TaskType>
bool PostProcessTask(const std::shared_ptr<TaskType> &task) {
  return task->PostProcessing();
}

template <class TaskType>
void ExpectValidation(const std::shared_ptr<TaskType> &task) {
  EXPECT_TRUE(ValidateTask(task));
}

template <class TaskType>
void ExpectPreProcessing(const std::shared_ptr<TaskType> &task) {
  EXPECT_TRUE(PreProcessTask(task));
}

template <class TaskType>
void ExpectRun(const std::shared_ptr<TaskType> &task) {
  EXPECT_TRUE(RunTask(task));
}

template <class TaskType>
void ExpectPostProcessing(const std::shared_ptr<TaskType> &task) {
  EXPECT_TRUE(PostProcessTask(task));
}

template <class TaskType>
void ExecuteTaskPipeline(const std::shared_ptr<TaskType> &task) {
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

template <class TaskType>
void RunAndCheckHull(const InType &points, size_t expected_size) {
  auto task = std::make_shared<TaskType>(points);
  ExecuteTaskPipeline(task);
  VerifyHull(task->GetOutput(), expected_size);
}

template <class TaskType>
void RunAndExpectFailure(const InType &points) {
  auto task = std::make_shared<TaskType>(points);
  EXPECT_FALSE(task->Validation());
  EXPECT_TRUE(task->GetOutput().empty());
}

void RunAndCheckAllImplementations(const InType &points, size_t expected_size) {
  RunAndCheckHull<UrinOGrahamPassageSEQ>(points, expected_size);
  RunAndCheckHull<UrinOGrahamPassageSTL>(points, expected_size);
  RunAndCheckHull<UrinOGrahamPassageOMP>(points, expected_size);
  RunAndCheckHull<UrinOGrahamPassageTBB>(points, expected_size);
  RunAndCheckHull<UrinOGrahamPassageALL>(points, expected_size);
}

void RunAndExpectAllImplementationsFailure(const InType &points) {
  RunAndExpectFailure<UrinOGrahamPassageSEQ>(points);
  RunAndExpectFailure<UrinOGrahamPassageSTL>(points);
  RunAndExpectFailure<UrinOGrahamPassageOMP>(points);
  RunAndExpectFailure<UrinOGrahamPassageTBB>(points);
  RunAndExpectFailure<UrinOGrahamPassageALL>(points);
}

TEST(UrinOGrahamPassage, EmptyInput) {
  RunAndExpectAllImplementationsFailure({});
}

TEST(UrinOGrahamPassage, SinglePoint) {
  RunAndExpectAllImplementationsFailure({Point(5.0, 3.0)});
}

TEST(UrinOGrahamPassage, TwoDistinctPoints) {
  RunAndExpectAllImplementationsFailure({Point(0.0, 0.0), Point(3.0, 4.0)});
}

TEST(UrinOGrahamPassage, CollinearPoints) {
  InType pts = {Point(0.0, 0.0), Point(1.0, 0.0), Point(2.0, 0.0), Point(3.0, 0.0), Point(4.0, 0.0)};
  RunAndCheckAllImplementations(pts, 2);
}

TEST(UrinOGrahamPassage, TrianglePoints) {
  InType pts = {Point(0.0, 0.0), Point(4.0, 0.0), Point(2.0, 3.0)};
  RunAndCheckAllImplementations(pts, 3);
}

TEST(UrinOGrahamPassage, SquarePoints) {
  InType pts = {Point(0.0, 0.0), Point(4.0, 0.0), Point(4.0, 4.0), Point(0.0, 4.0)};
  RunAndCheckAllImplementations(pts, 4);
}

TEST(UrinOGrahamPassage, SquareWithInteriorPoint) {
  InType pts = {Point(0.0, 0.0), Point(4.0, 0.0), Point(4.0, 4.0), Point(0.0, 4.0), Point(2.0, 2.0)};
  RunAndCheckAllImplementations(pts, 4);
}

TEST(UrinOGrahamPassage, RectangleWithCollinearPoints) {
  InType pts = {Point(0.0, 0.0), Point(1.0, 0.0), Point(2.0, 0.0), Point(3.0, 0.0),
                Point(3.0, 1.0), Point(2.0, 1.0), Point(1.0, 1.0), Point(0.0, 1.0)};
  RunAndCheckAllImplementations(pts, 4);
}

TEST(UrinOGrahamPassage, AllIdenticalPoints) {
  InType pts = {Point(3.0, 3.0), Point(3.0, 3.0), Point(3.0, 3.0), Point(3.0, 3.0), Point(3.0, 3.0)};
  RunAndExpectAllImplementationsFailure(pts);
}

TEST(UrinOGrahamPassage, PointOnBoundary) {
  InType pts = {Point(0.0, 0.0), Point(4.0, 0.0), Point(2.0, 0.0), Point(4.0, 4.0), Point(0.0, 4.0)};
  RunAndCheckAllImplementations(pts, 4);
}

TEST(UrinOGrahamPassage, VerticalCollinear) {
  InType pts = {Point(0.0, 0.0), Point(0.0, 1.0), Point(0.0, 2.0), Point(0.0, 5.0)};
  RunAndCheckAllImplementations(pts, 2);
}

TEST(UrinOGrahamPassage, LargeRandomSet) {
  InType pts;
  const int num_points = 100;
  pts.reserve(static_cast<size_t>(num_points));

  for (int i = 0; i < num_points; ++i) {
    const double angle = 2.0 * std::numbers::pi * static_cast<double>(i) / static_cast<double>(num_points);
    pts.emplace_back(std::cos(angle) * 10.0, std::sin(angle) * 10.0);
  }

  RunAndCheckAllImplementations(pts, 100);
}

TEST(UrinOGrahamPassage, HexagonWithCenter) {
  InType pts = {Point(2.0, 0.0),    Point(1.0, 1.73),  Point(-1.0, 1.73), Point(-2.0, 0.0),
                Point(-1.0, -1.73), Point(1.0, -1.73), Point(0.0, 0.0)};
  RunAndCheckAllImplementations(pts, 6);
}

}  // namespace
}  // namespace urin_o_graham_passage
