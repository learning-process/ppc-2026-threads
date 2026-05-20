#include "peterson_r_graham_scan_stl/stl/include/ops_stl.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <numbers>
#include <utility>
#include <vector>

#include "peterson_r_graham_scan_stl/common/include/common.hpp"

namespace peterson_r_graham_scan_stl {

namespace {
constexpr double kTolerance = 1e-12;

double CalculateOrientation(const Point2D &origin, const Point2D &a, const Point2D &b) {
  return ((a.coord_x - origin.coord_x) * (b.coord_y - origin.coord_y)) -
         ((a.coord_y - origin.coord_y) * (b.coord_x - origin.coord_x));
}

double CalculateSquaredDistance(const Point2D &first, const Point2D &second) {
  const double dx = first.coord_x - second.coord_x;
  const double dy = first.coord_y - second.coord_y;
  return (dx * dx) + (dy * dy);
}

class PointComparator {
 public:
  explicit PointComparator(const Point2D *reference) : origin_ptr_(reference) {}

  bool operator()(const Point2D &lhs, const Point2D &rhs) const {
    const double orientation = CalculateOrientation(*origin_ptr_, lhs, rhs);
    if (std::abs(orientation) > kTolerance) {
      return orientation > 0;
    }
    return CalculateSquaredDistance(*origin_ptr_, lhs) < CalculateSquaredDistance(*origin_ptr_, rhs);
  }

 private:
  const Point2D *origin_ptr_;
};
}  // namespace

PetersonGrahamScannerSTL::PetersonGrahamScannerSTL(const InputValue &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

void PetersonGrahamScannerSTL::LoadPoints(const PointSet &points) {
  input_points_ = points;
  external_data_provided_ = true;
}

PointSet PetersonGrahamScannerSTL::GetConvexHull() const {
  return hull_points_;
}

bool PetersonGrahamScannerSTL::ValidationImpl() {
  return GetInput() >= 0;
}

bool PetersonGrahamScannerSTL::PreProcessingImpl() {
  hull_points_.clear();

  if (!external_data_provided_) {
    input_points_.clear();
    const int count = GetInput();
    if (count <= 0) {
      return true;
    }

    input_points_.resize(count);
    const double angle_step = 2.0 * std::numbers::pi / count;

    for (int i = 0; i < count; ++i) {
      const double angle = angle_step * i;
      input_points_[i] = Point2D(std::cos(angle), std::sin(angle));
    }
  }

  return true;
}

bool PetersonGrahamScannerSTL::AreAllPointsIdentical(const PointSet &points) {
  if (points.empty()) {
    return true;
  }

  const Point2D &reference = points[0];
  return std::all_of(points.begin() + 1, points.end(), [&reference](const Point2D &point) {
    return std::abs(point.coord_x - reference.coord_x) <= kTolerance &&
           std::abs(point.coord_y - reference.coord_y) <= kTolerance;
  });
}

std::size_t PetersonGrahamScannerSTL::FindLowestPointParallel(const PointSet &points) {
  auto it = std::ranges::min_element(points, [](const Point2D &lhs, const Point2D &rhs) {
    if (lhs.coord_y < rhs.coord_y) {
      return true;
    }
    if (std::abs(lhs.coord_y - rhs.coord_y) < kTolerance) {
      return lhs.coord_x < rhs.coord_x;
    }
    return false;
  });

  return static_cast<std::size_t>(std::distance(points.begin(), it));
}

void PetersonGrahamScannerSTL::SortPointsByAngleParallel(PointSet &points) {
  if (points.size() < 2) {
    return;
  }

  const Point2D origin = points[0];
  PointComparator comparator(&origin);
  std::sort(points.begin() + 1, points.end(), comparator);
}

bool PetersonGrahamScannerSTL::RunImpl() {
  hull_points_.clear();
  const int total_points = static_cast<int>(input_points_.size());

  if (total_points == 0) {
    return true;
  }

  if (AreAllPointsIdentical(input_points_)) {
    hull_points_.push_back(input_points_.front());
    return true;
  }

  if (total_points < 3) {
    hull_points_ = input_points_;
    return true;
  }

  const std::size_t lowest_idx = FindLowestPointParallel(input_points_);
  std::swap(input_points_[0], input_points_[lowest_idx]);

  SortPointsByAngleParallel(input_points_);

  std::vector<Point2D> stack;
  stack.reserve(total_points);
  stack.push_back(input_points_[0]);
  stack.push_back(input_points_[1]);

  for (int i = 2; i < total_points; ++i) {
    while (static_cast<int>(stack.size()) >= 2) {
      const Point2D &second_last = stack[stack.size() - 2];
      const Point2D &last = stack.back();

      if (CalculateOrientation(second_last, last, input_points_[i]) <= kTolerance) {
        stack.pop_back();
      } else {
        break;
      }
    }
    stack.push_back(input_points_[i]);
  }

  hull_points_ = std::move(stack);
  return true;
}

bool PetersonGrahamScannerSTL::PostProcessingImpl() {
  GetOutput() = static_cast<int>(hull_points_.size());
  return true;
}

double PetersonGrahamScannerSTL::ComputeOrientation(const Point2D &origin, const Point2D &a, const Point2D &b) {
  return CalculateOrientation(origin, a, b);
}

double PetersonGrahamScannerSTL::ComputeDistanceSq(const Point2D &p1, const Point2D &p2) {
  return CalculateSquaredDistance(p1, p2);
}

}  // namespace peterson_r_graham_scan_stl
