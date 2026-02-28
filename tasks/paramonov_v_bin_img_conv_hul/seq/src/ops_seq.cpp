#include "paramonov_v_bin_img_conv_hull/seq/include/ops_seq.hpp"

#include <algorithm>
#include <deque>
#include <functional>
#include <iterator>
#include <limits>
#include <numeric>
#include <queue>
#include <set>
#include <stack>

namespace paramonov_v_bin_img_conv_hull {

namespace {
constexpr std::array<std::pair<int, int>, 4> kNeighbors = {{{1, 0}, {-1, 0}, {0, 1}, {0, -1}}};
}  // namespace

ConvexHullSequential::ConvexHullSequential(const InputType &input) {
  SetTypeOfTask(StaticTaskType());
  GetInput() = input;
}

bool ConvexHullSequential::ValidationImpl() {
  const auto &img = GetInput();
  if (img.rows <= 0 || img.cols <= 0) {
    return false;
  }

  const size_t expected_size = static_cast<size_t>(img.rows) * static_cast<size_t>(img.cols);
  return img.pixels.size() == expected_size;
}

bool ConvexHullSequential::PreProcessingImpl() {
  working_image_ = GetInput();
  BinarizeImage();
  GetOutput().clear();
  return true;
}

bool ConvexHullSequential::RunImpl() {
  ExtractConnectedComponents();
  return true;
}

bool ConvexHullSequential::PostProcessingImpl() {
  return true;
}

void ConvexHullSequential::BinarizeImage(uint8_t threshold) {
  std::transform(working_image_.pixels.begin(), working_image_.pixels.end(), working_image_.pixels.begin(),
                 [threshold](uint8_t pixel) { return pixel > threshold ? uint8_t{255} : uint8_t{0}; });
}

void ConvexHullSequential::FloodFill(int start_row, int start_col, std::vector<bool> &visited,
                                     std::vector<PixelPoint> &component) {
  std::stack<PixelPoint> pixel_stack;
  pixel_stack.emplace(start_row, start_col);

  const int rows = working_image_.rows;
  const int cols = working_image_.cols;

  visited[PixelIndex(start_row, start_col, cols)] = true;

  while (!pixel_stack.empty()) {
    PixelPoint current = pixel_stack.top();
    pixel_stack.pop();

    component.push_back(current);

    for (const auto &[dr, dc] : kNeighbors) {
      int nr = current.row + dr;
      int nc = current.col + dc;

      if (nr >= 0 && nr < rows && nc >= 0 && nc < cols) {
        size_t idx = PixelIndex(nr, nc, cols);
        if (!visited[idx] && working_image_.pixels[idx] == 255) {
          visited[idx] = true;
          pixel_stack.emplace(nr, nc);
        }
      }
    }
  }
}

void ConvexHullSequential::ExtractConnectedComponents() {
  const int rows = working_image_.rows;
  const int cols = working_image_.cols;
  const size_t total_pixels = static_cast<size_t>(rows) * static_cast<size_t>(cols);

  std::vector<bool> visited(total_pixels, false);
  auto &output_hulls = GetOutput();

  for (int r = 0; r < rows; ++r) {
    for (int c = 0; c < cols; ++c) {
      size_t idx = PixelIndex(r, c, cols);

      if (working_image_.pixels[idx] == 255 && !visited[idx]) {
        std::vector<PixelPoint> component;
        FloodFill(r, c, visited, component);

        if (!component.empty()) {
          std::vector<PixelPoint> hull = ComputeConvexHull(component);
          output_hulls.push_back(std::move(hull));
        }
      }
    }
  }
}

int64_t ConvexHullSequential::Orientation(const PixelPoint &p, const PixelPoint &q, const PixelPoint &r) {
  return (static_cast<int64_t>(q.col - p.col) * (r.row - p.row)) -
         (static_cast<int64_t>(q.row - p.row) * (r.col - p.col));
}

bool ConvexHullSequential::IsPointOnSegment(const PixelPoint &p, const PixelPoint &q, const PixelPoint &r) {
  if (Orientation(p, q, r) != 0) {
    return false;
  }

  return (q.col <= std::max(p.col, r.col) && q.col >= std::min(p.col, r.col) && q.row <= std::max(p.row, r.row) &&
          q.row >= std::min(p.row, r.row));
}

std::vector<PixelPoint> ConvexHullSequential::ComputeConvexHull(const std::vector<PixelPoint> &points) {
  if (points.size() <= 2) {
    return points;
  }

  auto lowest_point = *std::min_element(points.begin(), points.end(), [](const PixelPoint &a, const PixelPoint &b) {
    return (a.row == b.row) ? a.col < b.col : a.row < b.row;
  });

  std::vector<PixelPoint> sorted_points;
  std::copy_if(points.begin(), points.end(), std::back_inserter(sorted_points),
               [&lowest_point](const PixelPoint &p) { return !(p == lowest_point); });

  std::sort(sorted_points.begin(), sorted_points.end(), [&lowest_point](const PixelPoint &a, const PixelPoint &b) {
    int64_t orient = Orientation(lowest_point, a, b);
    if (orient == 0) {
      int64_t dist_a = (a.row - lowest_point.row) * (a.row - lowest_point.row) +
                       (a.col - lowest_point.col) * (a.col - lowest_point.col);
      int64_t dist_b = (b.row - lowest_point.row) * (b.row - lowest_point.row) +
                       (b.col - lowest_point.col) * (b.col - lowest_point.col);
      return dist_a < dist_b;
    }
    return orient > 0;
  });

  std::vector<PixelPoint> unique_points;
  unique_points.push_back(lowest_point);

  for (const auto &p : sorted_points) {
    while (unique_points.size() >= 2) {
      const auto &p1 = unique_points[unique_points.size() - 2];
      const auto &p2 = unique_points.back();

      if (Orientation(p1, p2, p) <= 0) {
        unique_points.pop_back();
      } else {
        break;
      }
    }
    unique_points.push_back(p);
  }

  std::vector<PixelPoint> hull;

  for (const auto &p : unique_points) {
    while (hull.size() >= 2) {
      const auto &a = hull[hull.size() - 2];
      const auto &b = hull.back();

      if (Orientation(a, b, p) <= 0) {
        hull.pop_back();
      } else {
        break;
      }
    }
    hull.push_back(p);
  }

  if (hull.size() > 1 && hull.front() == hull.back()) {
    hull.pop_back();
  }

  return hull;
}

}  // namespace paramonov_v_bin_img_conv_hull
