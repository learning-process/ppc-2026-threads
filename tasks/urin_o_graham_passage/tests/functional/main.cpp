#include <gtest/gtest.h>

#include <cmath>    // для std::cos, std::sin
#include <cstddef>  // для size_t
#include <memory>
#include <vector>

#include "urin_o_graham_passage/common/include/common.hpp"
#include "urin_o_graham_passage/seq/include/ops_seq.hpp"

namespace urin_o_graham_passage {
namespace {

bool IsConvexHull(const std::vector<Point> &hull) {
  if (hull.size() < 3) {
    return true;
  }

  for (size_t i = 0; i < hull.size(); ++i) {
    size_t prev = (i == 0) ? hull.size() - 1 : i - 1;
    size_t next = (i + 1) % hull.size();

    if (UrinOGrahamPassageSEQ::Orientation(hull[prev], hull[i], hull[next]) < 0) {
      return false;
    }
  }
  return true;
}

void RunAndCheckHull(const InType &points, size_t expected_size) {
  auto task = std::make_shared<UrinOGrahamPassageSEQ>(points);

  EXPECT_TRUE(task->Validation());
  EXPECT_TRUE(task->PreProcessing());
  EXPECT_TRUE(task->Run());
  EXPECT_TRUE(task->PostProcessing());

  const auto &hull = task->GetOutput();
  EXPECT_EQ(hull.size(), expected_size);
  EXPECT_TRUE(IsConvexHull(hull));
}

void RunAndExpectFailure(const InType &points) {
  auto task = std::make_shared<UrinOGrahamPassageSEQ>(points);
  EXPECT_FALSE(task->Validation());
  EXPECT_TRUE(task->GetOutput().empty());
}

TEST(UrinOGrahamPassageSeq, EmptyInput) {
  RunAndExpectFailure({});
}

TEST(UrinOGrahamPassageSeq, SinglePoint) {
  RunAndExpectFailure({Point(5.0, 3.0)});
}

TEST(UrinOGrahamPassageSeq, TwoDistinctPoints) {
  RunAndExpectFailure({Point(0.0, 0.0), Point(3.0, 4.0)});
}

TEST(UrinOGrahamPassageSeq, CollinearPoints) {
  InType pts = {Point(0.0, 0.0), Point(1.0, 0.0), Point(2.0, 0.0), Point(3.0, 0.0), Point(4.0, 0.0)};
  RunAndCheckHull(pts, 2);
}

TEST(UrinOGrahamPassageSeq, TrianglePoints) {
  InType pts = {Point(0.0, 0.0), Point(4.0, 0.0), Point(2.0, 3.0)};
  RunAndCheckHull(pts, 3);
}

TEST(UrinOGrahamPassageSeq, SquarePoints) {
  InType pts = {Point(0.0, 0.0), Point(4.0, 0.0), Point(4.0, 4.0), Point(0.0, 4.0)};
  RunAndCheckHull(pts, 4);
}

TEST(UrinOGrahamPassageSeq, SquareWithInteriorPoint) {
  InType pts = {Point(0.0, 0.0), Point(4.0, 0.0), Point(4.0, 4.0), Point(0.0, 4.0), Point(2.0, 2.0)};
  RunAndCheckHull(pts, 4);
}

TEST(UrinOGrahamPassageSeq, RectangleWithCollinearPoints) {
  InType pts = {Point(0.0, 0.0), Point(1.0, 0.0), Point(2.0, 0.0), Point(3.0, 0.0),
                Point(3.0, 1.0), Point(2.0, 1.0), Point(1.0, 1.0), Point(0.0, 1.0)};
  RunAndCheckHull(pts, 4);
}

TEST(UrinOGrahamPassageSeq, AllIdenticalPoints) {
  InType pts = {Point(3.0, 3.0), Point(3.0, 3.0), Point(3.0, 3.0), Point(3.0, 3.0), Point(3.0, 3.0)};
  RunAndExpectFailure(pts);
}

TEST(UrinOGrahamPassageSeq, PointOnBoundary) {
  InType pts = {Point(0.0, 0.0), Point(4.0, 0.0), Point(2.0, 0.0), Point(4.0, 4.0), Point(0.0, 4.0)};
  RunAndCheckHull(pts, 4);
}

TEST(UrinOGrahamPassageSeq, VerticalCollinear) {
  InType pts = {Point(0.0, 0.0), Point(0.0, 1.0), Point(0.0, 2.0), Point(0.0, 5.0)};
  RunAndCheckHull(pts, 2);
}

TEST(UrinOGrahamPassageSeq, LargeRandomSet) {
  InType pts;
  const int num_points = 100;
  pts.reserve(static_cast<size_t>(num_points));

  // const double pi = 3.14159265358979323846;  // Явно определяем pi
  for (int i = 0; i < num_points; ++i) {
    double angle = 2.0 * 3.14159 * static_cast<double>(i) / static_cast<double>(num_points);
    pts.emplace_back(std::cos(angle) * 10.0, std::sin(angle) * 10.0);
  }

  RunAndCheckHull(pts, 100);
}

TEST(UrinOGrahamPassageSeq, HexagonWithCenter) {
  InType pts = {Point(2.0, 0.0),    Point(1.0, 1.73),  Point(-1.0, 1.73), Point(-2.0, 0.0),
                Point(-1.0, -1.73), Point(1.0, -1.73), Point(0.0, 0.0)};
  RunAndCheckHull(pts, 6);
}

}  // namespace
}  // namespace urin_o_graham_passage
