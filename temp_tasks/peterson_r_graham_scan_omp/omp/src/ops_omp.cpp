#include "peterson_r_graham_scan_omp/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <utility>
#include <vector>

#include "peterson_r_graham_scan_omp/common/include/common.hpp"

namespace peterson_r_graham_scan_omp {

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

PetersonGrahamScannerOMP::PetersonGrahamScannerOMP(const InputValue &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

void PetersonGrahamScannerOMP::LoadPoints(const PointSet &points) {
  input_points_ = points;
  external_data_provided_ = true;
}

PointSet PetersonGrahamScannerOMP::GetConvexHull() const {
  return hull_points_;
}

bool PetersonGrahamScannerOMP::ValidationImpl() {
  return GetInput() >= 0;
}

bool PetersonGrahamScannerOMP::PreProcessingImpl() {
  hull_points_.clear();

  if (!external_data_provided_) {
    input_points_.clear();
    const int count = GetInput();
    if (count <= 0) {
      return true;
    }

    input_points_.resize(count);
    const double angle_step = 2.0 * std::numbers::pi / count;

#ifdef _MSC_VER
#  pragma omp parallel for
#else
#  pragma omp parallel for default(none) shared(input_points_, count, angle_step)
#endif
    for (int i = 0; i < count; ++i) {
      const double angle = angle_step * i;
      input_points_[i] = Point2D(std::cos(angle), std::sin(angle));
    }
  }

  return true;
}

bool PetersonGrahamScannerOMP::AreAllPointsIdentical(const PointSet &points) {
  if (points.empty()) {
    return true;
  }

  const Point2D &reference = points[0];
  const int num_points = static_cast<int>(points.size());

  bool all_identical = true;

#ifdef _MSC_VER
#  pragma omp parallel for reduction(&& : all_identical)
#else
#  pragma omp parallel for default(none) shared(points, reference, num_points, kTolerance) reduction(&& : all_identical)
#endif
  for (int i = 1; i < num_points; ++i) {
    if (std::abs(points[i].coord_x - reference.coord_x) > kTolerance ||
        std::abs(points[i].coord_y - reference.coord_y) > kTolerance) {
      all_identical = false;
    }
  }

  return all_identical;
}

std::size_t PetersonGrahamScannerOMP::FindLowestPointParallel(const PointSet &points) {
  const int total_points = static_cast<int>(points.size());
  int lowest_idx = 0;

#ifdef _MSC_VER
#  pragma omp parallel
#else
#  pragma omp parallel default(none) shared(points, total_points, kTolerance, lowest_idx)
#endif
  {
    int local_lowest = 0;

#pragma omp for nowait
    for (int i = 1; i < total_points; ++i) {
      if (points[i].coord_y < points[local_lowest].coord_y ||
          (std::abs(points[i].coord_y - points[local_lowest].coord_y) < kTolerance &&
           points[i].coord_x < points[local_lowest].coord_x)) {
        local_lowest = i;
      }
    }

#pragma omp critical
    {
      if (points[local_lowest].coord_y < points[lowest_idx].coord_y ||
          (std::abs(points[local_lowest].coord_y - points[lowest_idx].coord_y) < kTolerance &&
           points[local_lowest].coord_x < points[lowest_idx].coord_x)) {
        lowest_idx = local_lowest;
      }
    }
  }

  return static_cast<std::size_t>(lowest_idx);
}

void PetersonGrahamScannerOMP::ParallelMergeSort(PointSet &points, int left, int right, const Point2D &origin) {
  if (left >= right) {
    return;
  }

  int mid = left + ((right - left) / 2);

#ifdef _MSC_VER
#  pragma omp task shared(points)
#else
#  pragma omp task default(none) shared(points, left, mid, right, origin)
#endif
  ParallelMergeSort(points, left, mid, origin);

#ifdef _MSC_VER
#  pragma omp task shared(points)
#else
#  pragma omp task default(none) shared(points, left, mid, right, origin)
#endif
  ParallelMergeSort(points, mid + 1, right, origin);

#pragma omp taskwait

  PointComparator comp(&origin);
  std::inplace_merge(points.begin() + left, points.begin() + mid + 1, points.begin() + right + 1, comp);
}

void PetersonGrahamScannerOMP::SortPointsByAngleParallel(PointSet &points) {
  if (points.size() < 2) {
    return;
  }

  const Point2D origin = points[0];

#ifdef _MSC_VER
#  pragma omp parallel
#else
#  pragma omp parallel default(none) shared(points, origin)
#endif
  {
#pragma omp single
    {
      ParallelMergeSort(points, 1, static_cast<int>(points.size()) - 1, origin);
    }
  }
}

bool PetersonGrahamScannerOMP::RunImpl() {
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

bool PetersonGrahamScannerOMP::PostProcessingImpl() {
  GetOutput() = GetInput();
  return true;
}

double PetersonGrahamScannerOMP::ComputeOrientation(const Point2D &origin, const Point2D &a, const Point2D &b) {
  return CalculateOrientation(origin, a, b);
}

double PetersonGrahamScannerOMP::ComputeDistanceSq(const Point2D &p1, const Point2D &p2) {
  return CalculateSquaredDistance(p1, p2);
}

}  // namespace peterson_r_graham_scan_omp
