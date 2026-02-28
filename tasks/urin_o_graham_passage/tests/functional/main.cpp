#include <gtest/gtest.h>

#include <cstddef>
#include <memory>
#include <vector>

#include "urin_o_graham_passage/common/include/common.hpp"
#include "urin_o_graham_passage/seq/include/ops_seq.hpp"

namespace urin_o_graham_passage {
namespace {

static bool IsConvexHull(const std::vector<Point> &hull) {
  if (hull.size() < 3) {
    return true;
  }

  for (size_t i = 0; i < hull.size(); ++i) {
    const size_t prev = (i == 0) ? hull.size() - 1 : i - 1;
    const size_t next = (i + 1) % hull.size();

    if (UrinOGrahamPassageSEQ::Orientation(hull[prev], hull[i], hull[next]) < 0) {
      return false;
    }
  }
  return true;
}

TEST(UrinOGrahamPassageSeq, EmptyInput) {
  const InType empty_points;
  auto task = std::make_shared<UrinOGrahamPassageSEQ>(empty_points);

  EXPECT_FALSE(task->Validation());
  EXPECT_TRUE(task->GetOutput().empty());
}

TEST(UrinOGrahamPassageSeq, SinglePoint) {
  const InType pts = {Point(5.0, 3.0)};
  auto task = std::make_shared<UrinOGrahamPassageSEQ>(pts);

  EXPECT_FALSE(task->Validation());
  EXPECT_TRUE(task->GetOutput().empty());
}

TEST(UrinOGrahamPassageSeq, TwoDistinctPoints) {
  const InType pts = {Point(0.0, 0.0), Point(3.0, 4.0)};
  auto task = std::make_shared<UrinOGrahamPassageSEQ>(pts);

  EXPECT_FALSE(task->Validation());
  EXPECT_TRUE(task->GetOutput().empty());
}

TEST(UrinOGrahamPassageSeq, CollinearPoints) {
  const InType pts = {Point(0.0, 0.0), Point(1.0, 0.0), Point(2.0, 0.0), Point(3.0, 0.0), Point(4.0, 0.0)};
  auto task = std::make_shared<UrinOGrahamPassageSEQ>(pts);

  EXPECT_TRUE(task->Validation());
  EXPECT_TRUE(task->PreProcessing());
  EXPECT_TRUE(task->Run());
  EXPECT_TRUE(task->PostProcessing());

  const auto &hull = task->GetOutput();
  EXPECT_EQ(hull.size(), static_cast<size_t>(2));
  EXPECT_TRUE(IsConvexHull(hull));
}

TEST(UrinOGrahamPassageSeq, TrianglePoints) {
  const InType pts = {Point(0.0, 0.0), Point(4.0, 0.0), Point(2.0, 3.0)};
  auto task = std::make_shared<UrinOGrahamPassageSEQ>(pts);

  EXPECT_TRUE(task->Validation());
  EXPECT_TRUE(task->PreProcessing());
  EXPECT_TRUE(task->Run());
  EXPECT_TRUE(task->PostProcessing());

  const auto &hull = task->GetOutput();
  EXPECT_EQ(hull.size(), static_cast<size_t>(3));
  EXPECT_TRUE(IsConvexHull(hull));
}

TEST(UrinOGrahamPassageSeq, SquarePoints) {
  const InType pts = {Point(0.0, 0.0), Point(4.0, 0.0), Point(4.0, 4.0), Point(0.0, 4.0)};
  auto task = std::make_shared<UrinOGrahamPassageSEQ>(pts);

  EXPECT_TRUE(task->Validation());
  EXPECT_TRUE(task->PreProcessing());
  EXPECT_TRUE(task->Run());
  EXPECT_TRUE(task->PostProcessing());

  const auto &hull = task->GetOutput();
  EXPECT_EQ(hull.size(), static_cast<size_t>(4));
  EXPECT_TRUE(IsConvexHull(hull));
}

TEST(UrinOGrahamPassageSeq, SquareWithInteriorPoint) {
  const InType pts = {Point(0.0, 0.0), Point(4.0, 0.0), Point(4.0, 4.0), Point(0.0, 4.0), Point(2.0, 2.0)};
  auto task = std::make_shared<UrinOGrahamPassageSEQ>(pts);

  EXPECT_TRUE(task->Validation());
  EXPECT_TRUE(task->PreProcessing());
  EXPECT_TRUE(task->Run());
  EXPECT_TRUE(task->PostProcessing());

  const auto &hull = task->GetOutput();
  EXPECT_EQ(hull.size(), static_cast<size_t>(4));
  EXPECT_TRUE(IsConvexHull(hull));
}

TEST(UrinOGrahamPassageSeq, RectangleWithCollinearPoints) {
  const InType pts = {Point(0.0, 0.0), Point(1.0, 0.0), Point(2.0, 0.0), Point(3.0, 0.0),
                      Point(3.0, 1.0), Point(2.0, 1.0), Point(1.0, 1.0), Point(0.0, 1.0)};
  auto task = std::make_shared<UrinOGrahamPassageSEQ>(pts);

  EXPECT_TRUE(task->Validation());
  EXPECT_TRUE(task->PreProcessing());
  EXPECT_TRUE(task->Run());
  EXPECT_TRUE(task->PostProcessing());

  const auto &hull = task->GetOutput();
  EXPECT_EQ(hull.size(), static_cast<size_t>(4));
  EXPECT_TRUE(IsConvexHull(hull));
}

TEST(UrinOGrahamPassageSeq, AllIdenticalPoints) {
  const InType pts = {Point(3.0, 3.0), Point(3.0, 3.0), Point(3.0, 3.0), Point(3.0, 3.0), Point(3.0, 3.0)};
  auto task = std::make_shared<UrinOGrahamPassageSEQ>(pts);

  EXPECT_FALSE(task->Validation());
  EXPECT_TRUE(task->GetOutput().empty());
}

TEST(UrinOGrahamPassageSeq, PointOnBoundary) {
  const InType pts = {Point(0.0, 0.0), Point(4.0, 0.0), Point(2.0, 0.0), Point(4.0, 4.0), Point(0.0, 4.0)};
  auto task = std::make_shared<UrinOGrahamPassageSEQ>(pts);

  EXPECT_TRUE(task->Validation());
  EXPECT_TRUE(task->PreProcessing());
  EXPECT_TRUE(task->Run());
  EXPECT_TRUE(task->PostProcessing());

  const auto &hull = task->GetOutput();
  EXPECT_EQ(hull.size(), static_cast<size_t>(4));
  EXPECT_TRUE(IsConvexHull(hull));
}

TEST(UrinOGrahamPassageSeq, VerticalCollinear) {
  const InType pts = {Point(0.0, 0.0), Point(0.0, 1.0), Point(0.0, 2.0), Point(0.0, 5.0)};
  auto task = std::make_shared<UrinOGrahamPassageSEQ>(pts);

  EXPECT_TRUE(task->Validation());
  EXPECT_TRUE(task->PreProcessing());
  EXPECT_TRUE(task->Run());
  EXPECT_TRUE(task->PostProcessing());

  const auto &hull = task->GetOutput();
  EXPECT_EQ(hull.size(), static_cast<size_t>(2));
  EXPECT_TRUE(IsConvexHull(hull));
}

TEST(UrinOGrahamPassageSeq, LargeRandomSet) {
  InType pts;
  const int num_points = 100;
  pts.reserve(static_cast<size_t>(num_points));

  for (int i = 0; i < num_points; ++i) {
    const double angle = 2.0 * 3.14159 * static_cast<double>(i) / static_cast<double>(num_points);
    pts.emplace_back(std::cos(angle) * 10.0, std::sin(angle) * 10.0);
  }

  auto task = std::make_shared<UrinOGrahamPassageSEQ>(pts);

  EXPECT_TRUE(task->Validation());
  EXPECT_TRUE(task->PreProcessing());
  EXPECT_TRUE(task->Run());
  EXPECT_TRUE(task->PostProcessing());

  const auto &hull = task->GetOutput();
  EXPECT_GE(hull.size(), static_cast<size_t>(3));
  EXPECT_TRUE(IsConvexHull(hull));
}

TEST(UrinOGrahamPassageSeq, HexagonWithCenter) {
  const InType pts = {Point(2.0, 0.0),    Point(1.0, 1.73),  Point(-1.0, 1.73), Point(-2.0, 0.0),
                      Point(-1.0, -1.73), Point(1.0, -1.73), Point(0.0, 0.0)};
  auto task = std::make_shared<UrinOGrahamPassageSEQ>(pts);

  EXPECT_TRUE(task->Validation());
  EXPECT_TRUE(task->PreProcessing());
  EXPECT_TRUE(task->Run());
  EXPECT_TRUE(task->PostProcessing());

  const auto &hull = task->GetOutput();
  EXPECT_EQ(hull.size(), static_cast<size_t>(6));
  EXPECT_TRUE(IsConvexHull(hull));
}

}  // namespace
}  // namespace urin_o_graham_passage
