#include "shkenev_i_constr_hull_for_binary_image_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <queue>
#include <vector>

namespace shkenev_i_constr_hull_for_binary_image_seq {

namespace {

const std::array<std::pair<int, int>, 4> kNeighborDirections = {std::make_pair(1, 0), std::make_pair(-1, 0),
                                                                std::make_pair(0, 1), std::make_pair(0, -1)};

void ExploreRegion(int start_row, int start_col, const BinaryImageData &image, std::vector<bool> &visited_mask,
                   std::vector<ImagePoint> &region) {
  std::queue<ImagePoint> exploration_queue;
  exploration_queue.emplace(start_row, start_col);
  visited_mask[helpers::ComputeIndex(start_row, start_col, image.cols)] = true;

  while (!exploration_queue.empty()) {
    ImagePoint current = exploration_queue.front();
    exploration_queue.pop();
    region.push_back(current);

    for (const auto &[dr, dc] : kNeighborDirections) {
      int neighbor_row = current.row + dr;
      int neighbor_col = current.col + dc;

      if (!helpers::IsValidCoordinate(neighbor_row, neighbor_col, image.rows, image.cols)) {
        continue;
      }

      size_t neighbor_idx = helpers::ComputeIndex(neighbor_row, neighbor_col, image.cols);
      if (visited_mask[neighbor_idx] || image.data[neighbor_idx] == 0) {
        continue;
      }

      visited_mask[neighbor_idx] = true;
      exploration_queue.emplace(neighbor_row, neighbor_col);
    }
  }
}

}  // namespace

ShkenevIConstrHullSeq::ShkenevIConstrHullSeq(const InType &in) : working_image_(in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool ShkenevIConstrHullSeq::ValidationImpl() {
  const auto &input = GetInput();

  if (input.rows <= 0 || input.cols <= 0) {
    return false;
  }

  size_t expected_size = static_cast<size_t>(input.rows) * static_cast<size_t>(input.cols);
  return input.data.size() == expected_size;
}

bool ShkenevIConstrHullSeq::PreProcessingImpl() {
  working_image_ = GetInput();
  BinarizeImage();
  return true;
}

bool ShkenevIConstrHullSeq::RunImpl() {
  ExtractConnectedRegions();

  working_image_.convex_envelopes.clear();
  working_image_.convex_envelopes.reserve(working_image_.connected_components.size());

  for (const auto &region : working_image_.connected_components) {
    if (region.empty()) {
      continue;
    }

    if (region.size() <= 2) {
      working_image_.convex_envelopes.push_back(region);
    } else {
      working_image_.convex_envelopes.push_back(ComputeConvexEnvelope(region));
    }
  }

  GetOutput() = working_image_;
  return true;
}

bool ShkenevIConstrHullSeq::PostProcessingImpl() {
  return true;
}

void ShkenevIConstrHullSeq::BinarizeImage(uint8_t threshold) {
  for (auto &pixel : working_image_.data) {
    pixel = (pixel > threshold) ? static_cast<uint8_t>(1) : static_cast<uint8_t>(0);
  }
}

void ShkenevIConstrHullSeq::ExtractConnectedRegions() {
  size_t total_pixels = static_cast<size_t>(working_image_.rows) * static_cast<size_t>(working_image_.cols);
  std::vector<bool> visited(total_pixels, false);

  working_image_.connected_components.clear();

  for (int r = 0; r < working_image_.rows; ++r) {
    for (int c = 0; c < working_image_.cols; ++c) {
      size_t idx = helpers::ComputeIndex(r, c, working_image_.cols);

      if (working_image_.data[idx] == 0 || visited[idx]) {
        continue;
      }

      std::vector<ImagePoint> region;
      ExploreRegion(r, c, working_image_, visited, region);

      if (!region.empty()) {
        working_image_.connected_components.push_back(std::move(region));
      }
    }
  }
}

long long ShkenevIConstrHullSeq::ComputeCrossProduct(const ImagePoint &p1, const ImagePoint &p2, const ImagePoint &p3) {
  long long v1_row = static_cast<long long>(p2.row) - static_cast<long long>(p1.row);
  long long v1_col = static_cast<long long>(p2.col) - static_cast<long long>(p1.col);
  long long v2_row = static_cast<long long>(p3.row) - static_cast<long long>(p2.row);
  long long v2_col = static_cast<long long>(p3.col) - static_cast<long long>(p2.col);

  return (v1_row * v2_col) - (v1_col * v2_row);
}

std::vector<ImagePoint> ShkenevIConstrHullSeq::ComputeConvexEnvelope(const std::vector<ImagePoint> &region) {
  if (region.size() <= 2) {
    return region;
  }

  std::vector<ImagePoint> sorted_points = region;

  std::sort(sorted_points.begin(), sorted_points.end(), [](const ImagePoint &a, const ImagePoint &b) {
    return (a.row != b.row) ? (a.row < b.row) : (a.col < b.col);
  });

  auto last_unique = std::unique(sorted_points.begin(), sorted_points.end());
  sorted_points.erase(last_unique, sorted_points.end());

  if (sorted_points.size() <= 2) {
    return sorted_points;
  }

  std::deque<ImagePoint> lower_hull;
  std::deque<ImagePoint> upper_hull;

  for (const auto &point : sorted_points) {
    while (lower_hull.size() >= 2) {
      const auto &p1 = lower_hull[lower_hull.size() - 2];
      const auto &p2 = lower_hull.back();

      if (ComputeCrossProduct(p1, p2, point) <= 0) {
        lower_hull.pop_back();
      } else {
        break;
      }
    }
    lower_hull.push_back(point);
  }

  for (auto it = sorted_points.rbegin(); it != sorted_points.rend(); ++it) {
    const auto &point = *it;

    while (upper_hull.size() >= 2) {
      const auto &p1 = upper_hull[upper_hull.size() - 2];
      const auto &p2 = upper_hull.back();

      if (ComputeCrossProduct(p1, p2, point) <= 0) {
        upper_hull.pop_back();
      } else {
        break;
      }
    }
    upper_hull.push_back(point);
  }

  lower_hull.pop_back();
  upper_hull.pop_back();

  std::vector<ImagePoint> convex_envelope;
  convex_envelope.reserve(lower_hull.size() + upper_hull.size());

  convex_envelope.insert(convex_envelope.end(), lower_hull.begin(), lower_hull.end());
  convex_envelope.insert(convex_envelope.end(), upper_hull.begin(), upper_hull.end());

  return convex_envelope;
}

}  // namespace shkenev_i_constr_hull_for_binary_image_seq
