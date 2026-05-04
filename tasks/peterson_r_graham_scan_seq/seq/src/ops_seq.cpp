#include "peterson_r_graham_scan_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <utility>
#include <vector>

#include "peterson_r_graham_scan_seq/common/include/common.hpp"

namespace peterson_r_graham_scan_seq {

namespace {
constexpr double kTolerance = 1e-12;

double CalculateOrientation(const Point2D &origin, const Point2D &a, const Point2D &b) {
  // Добавлены скобки для явного указания порядка операций
  return ((a.coord_x - origin.coord_x) * (b.coord_y - origin.coord_y)) -
         ((a.coord_y - origin.coord_y) * (b.coord_x - origin.coord_x));
}

double CalculateSquaredDistance(const Point2D &first, const Point2D &second) {
  const double dx = first.coord_x - second.coord_x;
  const double dy = first.coord_y - second.coord_y;
  // Добавлены скобки для явного указания порядка операций
  return (dx * dx) + (dy * dy);
}

// Убрана ссылка на const в члене класса
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
  const Point2D *origin_ptr_;  // Используем указатель вместо ссылки
};
}  // namespace

PetersonGrahamScanner::PetersonGrahamScanner(const InputValue &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

void PetersonGrahamScanner::LoadPoints(const PointSet &points) {
  input_points_ = points;
  external_data_provided_ = true;
}

PointSet PetersonGrahamScanner::GetConvexHull() const {
  return hull_points_;
}

bool PetersonGrahamScanner::ValidationImpl() {
  return GetInput() >= 0;
}

bool PetersonGrahamScanner::PreProcessingImpl() {
  hull_points_.clear();

  if (!external_data_provided_) {
    input_points_.clear();
    const int count = GetInput();
    if (count <= 0) {
      return true;
    }

    input_points_.reserve(count);
    const double angle_step = 2.0 * std::numbers::pi / count;

    for (int i = 0; i < count; ++i) {
      const double angle = angle_step * i;
      input_points_.emplace_back(std::cos(angle), std::sin(angle));
    }
  }

  return true;
}

bool PetersonGrahamScanner::RunImpl() {
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

  // Find the point with minimum y (and leftmost if tie)
  const std::size_t lowest_idx = FindLowestPoint(input_points_);
  std::swap(input_points_[0], input_points_[lowest_idx]);

  // Sort remaining points by polar angle
  SortByAngle(input_points_);

  // Graham scan algorithm
  std::vector<Point2D> stack;
  stack.reserve(total_points);
  stack.push_back(input_points_[0]);
  stack.push_back(input_points_[1]);

  for (int i = 2; i < total_points; ++i) {
    while (stack.size() >= 2) {
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

bool PetersonGrahamScanner::PostProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

double PetersonGrahamScanner::ComputeOrientation(const Point2D &origin, const Point2D &a, const Point2D &b) {
  return CalculateOrientation(origin, a, b);
}

double PetersonGrahamScanner::ComputeDistanceSq(const Point2D &p1, const Point2D &p2) {
  return CalculateSquaredDistance(p1, p2);
}

bool PetersonGrahamScanner::AreAllPointsIdentical(const PointSet &points) {
  if (points.empty()) {
    return true;
  }

  const Point2D &reference = points[0];

  for (std::size_t i = 1; i < points.size(); ++i) {
    if (std::abs(points[i].coord_x - reference.coord_x) > kTolerance ||
        std::abs(points[i].coord_y - reference.coord_y) > kTolerance) {
      return false;
    }
  }

  return true;
}

std::size_t PetersonGrahamScanner::FindLowestPoint(const PointSet &points) {
  std::size_t lowest = 0;

  for (std::size_t i = 1; i < points.size(); ++i) {
    if (points[i].coord_y < points[lowest].coord_y ||
        (std::abs(points[i].coord_y - points[lowest].coord_y) < kTolerance &&
         points[i].coord_x < points[lowest].coord_x)) {
      lowest = i;
    }
  }

  return lowest;
}

void PetersonGrahamScanner::SortByAngle(PointSet &points) {
  if (points.size() < 2) {
    return;
  }

  const Point2D origin = points[0];
  PointComparator comparator(&origin);  // Передаем указатель
  std::sort(points.begin() + 1, points.end(), comparator);
}

}  // namespace peterson_r_graham_scan_seq
