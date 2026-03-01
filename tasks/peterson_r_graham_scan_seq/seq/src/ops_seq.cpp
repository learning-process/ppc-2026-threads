#include "peterson_r_graham_scan_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>

#include "peterson_r_graham_scan_seq/common/include/common.hpp"

namespace peterson_r_graham_scan_seq {

namespace {
constexpr double kEpsilon = 1e-12;
}  // namespace

PetersonRGrahamScanSeq::PetersonRGrahamScanSeq(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = 0;
}

void PetersonRGrahamScanSeq::SetPoints(const PointCloud &cloud) {
  dataset_ = cloud;
  user_provided_data_ = true;
}

PointCloud PetersonRGrahamScanSeq::GetHull() const {
  return boundary_;
}

bool PetersonRGrahamScanSeq::ValidationImpl() {
  return GetInput() >= 0;
}

bool PetersonRGrahamScanSeq::PreProcessingImpl() {
  boundary_.clear();

  if (!user_provided_data_) {
    dataset_.clear();
    if (GetInput() <= 0) {
      return true;
    }

    const int sample_size = GetInput();
    dataset_.reserve(sample_size);

    for (int idx = 0; idx < sample_size; ++idx) {
      const double phase = (2.0 * std::numbers::pi * static_cast<double>(idx)) / static_cast<double>(sample_size);
      dataset_.emplace_back(std::cos(phase), std::sin(phase));
    }
  } else {
    if (static_cast<std::size_t>(GetInput()) != dataset_.size()) {
      GetInput() = static_cast<InType>(dataset_.size());
    }
  }

  return true;
}

bool PetersonRGrahamScanSeq::RunImpl() {
  boundary_.clear();
  const int total = static_cast<int>(dataset_.size());

  if (total == 0) {
    return true;
  }

  if (IsUniformCloud(dataset_)) {
    boundary_.push_back(dataset_.front());
    return true;
  }

  if (total < 3) {
    boundary_ = dataset_;
    return true;
  }

  const std::size_t anchor_pos = FindAnchorIndex(dataset_);
  std::iter_swap(dataset_.begin(), dataset_.begin() + static_cast<std::ptrdiff_t>(anchor_pos));

  ArrangeByPolarAngle(dataset_);

  boundary_.reserve(total);
  boundary_.push_back(dataset_[0]);
  boundary_.push_back(dataset_[1]);

  for (int i = 2; i < total; ++i) {
    while (boundary_.size() >= 2) {
      const Coordinate2D &prev2 = boundary_[boundary_.size() - 2];
      const Coordinate2D &prev1 = boundary_.back();

      if (CrossProduct(prev2, prev1, dataset_[i]) <= kEpsilon) {
        boundary_.pop_back();
      } else {
        break;
      }
    }
    boundary_.push_back(dataset_[i]);
  }

  return true;
}

bool PetersonRGrahamScanSeq::PostProcessingImpl() {
  GetOutput() = static_cast<OutType>(boundary_.size());
  return true;
}

double PetersonRGrahamScanSeq::CrossProduct(const Coordinate2D &origin, const Coordinate2D &a, const Coordinate2D &b) {
  return ((a.x - origin.x) * (b.y - origin.y)) - ((a.y - origin.y) * (b.x - origin.x));
}

double PetersonRGrahamScanSeq::SquaredDistance(const Coordinate2D &p1, const Coordinate2D &p2) {
  const double dx = p1.x - p2.x;
  const double dy = p1.y - p2.y;
  return (dx * dx) + (dy * dy);
}

bool PetersonRGrahamScanSeq::IsUniformCloud(const PointCloud &cloud) {
  if (cloud.empty()) {
    return true;
  }

  const Coordinate2D reference = cloud[0];

  for (std::size_t i = 1; i < cloud.size(); ++i) {
    if (std::abs(cloud[i].x - reference.x) > kEpsilon || std::abs(cloud[i].y - reference.y) > kEpsilon) {
      return false;
    }
  }

  return true;
}

std::size_t PetersonRGrahamScanSeq::FindAnchorIndex(const PointCloud &cloud) {
  std::size_t min_idx = 0;

  for (std::size_t i = 1; i < cloud.size(); ++i) {
    if (cloud[i].y < cloud[min_idx].y ||
        (std::abs(cloud[i].y - cloud[min_idx].y) < kEpsilon && cloud[i].x < cloud[min_idx].x)) {
      min_idx = i;
    }
  }

  return min_idx;
}

void PetersonRGrahamScanSeq::ArrangeByPolarAngle(PointCloud &cloud) {
  if (cloud.size() < 2) {
    return;
  }

  const Coordinate2D origin = cloud[0];

  auto comparator = [&origin](const Coordinate2D &lhs, const Coordinate2D &rhs) -> bool {
    const double orient = CrossProduct(origin, lhs, rhs);
    if (orient > kEpsilon) {
      return true;
    }
    if (orient < -kEpsilon) {
      return false;
    }
    return SquaredDistance(origin, lhs) < SquaredDistance(origin, rhs);
  };

  std::sort(cloud.begin() + 1, cloud.end(), comparator);
}

}  // namespace peterson_r_graham_scan_seq
