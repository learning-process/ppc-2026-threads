#include "peterson_r_graham_scan_tbb/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <atomic>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <utility>
#include <vector>

#include "peterson_r_graham_scan_tbb/common/include/common.hpp"

namespace peterson_r_graham_scan_tbb {

namespace {
constexpr double kTolerance = 0.0;

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

PetersonGrahamScannerTBB::PetersonGrahamScannerTBB(const InputValue &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

void PetersonGrahamScannerTBB::LoadPoints(const PointSet &points) {
  input_points_ = points;
  external_data_provided_ = true;
}

PointSet PetersonGrahamScannerTBB::GetConvexHull() const {
  return hull_points_;
}

bool PetersonGrahamScannerTBB::ValidationImpl() {
  return GetInput() >= 0;
}

bool PetersonGrahamScannerTBB::PreProcessingImpl() {
  hull_points_.clear();

  if (!external_data_provided_) {
    input_points_.clear();
    const int count = GetInput();
    if (count <= 0) {
      return true;
    }

    input_points_.resize(count);
    const double angle_step = 2.0 * std::numbers::pi / count;

    tbb::parallel_for(tbb::blocked_range<int>(0, count), [this, angle_step](const tbb::blocked_range<int> &range) {
      for (int i = range.begin(); i != range.end(); ++i) {
        const double angle = angle_step * i;
        input_points_[i] = Point2D(std::cos(angle), std::sin(angle));
      }
    });
  }

  return true;
}

bool PetersonGrahamScannerTBB::AreAllPointsIdentical(const PointSet &points) {
  if (points.empty()) {
    return true;
  }

  const Point2D &reference = points[0];
  std::atomic<bool> all_identical{true};

  tbb::parallel_for(tbb::blocked_range<std::size_t>(1, points.size()),
                    [&points, &reference, &all_identical](const tbb::blocked_range<std::size_t> &range) {
    for (std::size_t i = range.begin(); i != range.end(); ++i) {
      if (std::abs(points[i].coord_x - reference.coord_x) > kTolerance ||
          std::abs(points[i].coord_y - reference.coord_y) > kTolerance) {
        all_identical.store(false);
      }
    }
  });

  return all_identical.load();
}

std::size_t PetersonGrahamScannerTBB::FindLowestPointParallel(const PointSet &points) {
  struct LowestPointResult {
    std::size_t index = 0;
    double min_y = 0.0;
    double min_x = 0.0;

    LowestPointResult() = default;
    LowestPointResult(std::size_t idx, double y, double x) : index(idx), min_y(y), min_x(x) {}
  };

  LowestPointResult result = tbb::parallel_reduce(
      tbb::blocked_range<std::size_t>(0, points.size()), LowestPointResult(0, points[0].coord_y, points[0].coord_x),
      [&points](const tbb::blocked_range<std::size_t> &range, LowestPointResult current_result) {
    for (std::size_t i = range.begin(); i != range.end(); ++i) {
      if (points[i].coord_y < current_result.min_y ||
          (std::abs(points[i].coord_y - current_result.min_y) < kTolerance &&
           points[i].coord_x < current_result.min_x)) {
        current_result = LowestPointResult(i, points[i].coord_y, points[i].coord_x);
      }
    }
    return current_result;
  }, [](LowestPointResult lhs, LowestPointResult rhs) {
    if (rhs.min_y < lhs.min_y || (std::abs(rhs.min_y - lhs.min_y) < kTolerance && rhs.min_x < lhs.min_x)) {
      return rhs;
    }
    return lhs;
  });

  return result.index;
}

void PetersonGrahamScannerTBB::SortPointsByAngleParallel(PointSet &points) {
  if (points.size() < 2) {
    return;
  }

  const Point2D origin = points[0];
  PointComparator comparator(&origin);
  tbb::parallel_sort(points.begin() + 1, points.end(), comparator);
}

bool PetersonGrahamScannerTBB::RunImpl() {
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

bool PetersonGrahamScannerTBB::PostProcessingImpl() {
  GetOutput() = static_cast<int>(hull_points_.size());
  return true;
}

double PetersonGrahamScannerTBB::ComputeOrientation(const Point2D &origin, const Point2D &a, const Point2D &b) {
  return CalculateOrientation(origin, a, b);
}

double PetersonGrahamScannerTBB::ComputeDistanceSq(const Point2D &p1, const Point2D &p2) {
  return CalculateSquaredDistance(p1, p2);
}

}  // namespace peterson_r_graham_scan_tbb
