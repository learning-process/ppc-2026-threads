#include "paramonov_v_bin_img_conv_hul_stl/stl/include/ops_stl.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <stack>
#include <utility>
#include <vector>

#include "paramonov_v_bin_img_conv_hul_stl/common/include/common.hpp"

namespace paramonov_v_bin_img_conv_hul_stl {

namespace {
constexpr std::array<std::pair<int, int>, 4> kNeighbors = {{{1, 0}, {-1, 0}, {0, 1}, {0, -1}}};

bool ComparePoints(const PixelPoint &a, const PixelPoint &b) {
  if (a.row != b.row) {
    return a.row < b.row;
  }
  return a.col < b.col;
}

int64_t SquaredDistance(const PixelPoint &a, const PixelPoint &b) {
  auto dx = static_cast<int64_t>(a.col - b.col);
  auto dy = static_cast<int64_t>(a.row - b.row);
  return (dx * dx) + (dy * dy);
}

}  // namespace

ConvexHullSTL::ConvexHullSTL(const InputType &input) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = input;
}

bool ConvexHullSTL::ValidationImpl() {
  const auto &img = GetInput();
  if (img.rows <= 0 || img.cols <= 0) {
    return false;
  }

  const size_t expected_size = static_cast<size_t>(img.rows) * static_cast<size_t>(img.cols);
  return img.pixels.size() == expected_size;
}

bool ConvexHullSTL::PreProcessingImpl() {
  working_image_ = GetInput();
  BinarizeImage();
  GetOutput().clear();
  return true;
}

bool ConvexHullSTL::RunImpl() {
  ExtractConnectedComponents();
  return true;
}

bool ConvexHullSTL::PostProcessingImpl() {
  return true;
}

void ConvexHullSTL::BinarizeImage(uint8_t threshold) {
  const size_t size = working_image_.pixels.size();
  auto &pixels = working_image_.pixels;

  for (size_t i = 0; i < size; ++i) {
    pixels[i] = pixels[i] > threshold ? uint8_t{255} : uint8_t{0};
  }
}

void ConvexHullSTL::FloodFill(int start_row, int start_col, std::vector<bool> &visited,
                              std::vector<PixelPoint> &component) const {
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
      int next_row = current.row + dr;
      int next_col = current.col + dc;

      if (next_row >= 0 && next_row < rows && next_col >= 0 && next_col < cols) {
        size_t idx = PixelIndex(next_row, next_col, cols);
        if (!visited[idx] && working_image_.pixels[idx] == 255) {
          visited[idx] = true;
          pixel_stack.emplace(next_row, next_col);
        }
      }
    }
  }
}

void ConvexHullSTL::ExtractConnectedComponents() {
  const int rows = working_image_.rows;
  const int cols = working_image_.cols;
  const size_t total_pixels = static_cast<size_t>(rows) * static_cast<size_t>(cols);

  std::vector<bool> visited(total_pixels, false);
  auto &output_hulls = GetOutput();

  for (int row = 0; row < rows; ++row) {
    for (int col = 0; col < cols; ++col) {
      size_t idx = PixelIndex(row, col, cols);

      if (working_image_.pixels[idx] == 255 && !visited[idx]) {
        std::vector<PixelPoint> component;
        FloodFill(row, col, visited, component);

        if (!component.empty()) {
          std::vector<PixelPoint> hull = ComputeConvexHull(component);
          output_hulls.push_back(std::move(hull));
        }
      }
    }
  }
}

int64_t ConvexHullSTL::Orientation(const PixelPoint &p, const PixelPoint &q, const PixelPoint &r) {
  return (static_cast<int64_t>(q.col - p.col) * (r.row - p.row)) -
         (static_cast<int64_t>(q.row - p.row) * (r.col - p.col));
}

PixelPoint ConvexHullSTL::FindLowestPoint(const std::vector<PixelPoint> &points) {
  return *std::ranges::min_element(points, ComparePoints);
}

std::vector<PixelPoint> ConvexHullSTL::SortPointsByAngle(const std::vector<PixelPoint> &points,
                                                         const PixelPoint &lowest_point) {
  std::vector<PixelPoint> sorted_points;
  sorted_points.reserve(points.size() - 1);

  std::ranges::copy_if(points, std::back_inserter(sorted_points), [&lowest_point](const PixelPoint &p) {
    return (p.row != lowest_point.row) || (p.col != lowest_point.col);
  });

  std::ranges::sort(sorted_points, [&lowest_point](const PixelPoint &a, const PixelPoint &b) {
    int64_t orient = Orientation(lowest_point, a, b);
    if (orient != 0) {
      return orient > 0;
    }
    return SquaredDistance(a, lowest_point) < SquaredDistance(b, lowest_point);
  });

  return sorted_points;
}

std::vector<PixelPoint> ConvexHullSTL::RemoveCollinearPoints(const std::vector<PixelPoint> &sorted_points,
                                                             const PixelPoint &lowest_point) {
  if (sorted_points.empty()) {
    return {};
  }

  std::vector<PixelPoint> unique_points;
  unique_points.reserve(sorted_points.size());

  for (size_t i = 0; i < sorted_points.size(); ++i) {
    bool is_last = (i == sorted_points.size() - 1);
    bool is_collinear = !is_last && Orientation(lowest_point, sorted_points[i], sorted_points[i + 1]) == 0;

    if (is_last || !is_collinear) {
      unique_points.push_back(sorted_points[i]);
    }
  }

  return unique_points;
}

std::vector<PixelPoint> ConvexHullSTL::HandleCollinearCase(const std::vector<PixelPoint> &points,
                                                           const PixelPoint &lowest_point) {
  std::vector<PixelPoint> collinear_hull;
  collinear_hull.push_back(lowest_point);

  if (!points.empty()) {
    collinear_hull.push_back(points.back());
  }

  if (collinear_hull.size() == 2) {
    bool need_swap = (collinear_hull[0].row > collinear_hull[1].row) ||
                     (collinear_hull[0].row == collinear_hull[1].row && collinear_hull[0].col > collinear_hull[1].col);
    if (need_swap) {
      std::swap(collinear_hull[0], collinear_hull[1]);
    }
  }

  return collinear_hull;
}

std::vector<PixelPoint> ConvexHullSTL::BuildHull(const std::vector<PixelPoint> &unique_points,
                                                 const PixelPoint &lowest_point) {
  std::vector<PixelPoint> hull;
  hull.reserve(unique_points.size() + 1);
  hull.push_back(lowest_point);
  hull.push_back(unique_points[0]);

  for (size_t i = 1; i < unique_points.size(); ++i) {
    const auto &p = unique_points[i];
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

  return hull;
}

std::vector<PixelPoint> ConvexHullSTL::ComputeConvexHull(const std::vector<PixelPoint> &points) {
  if (points.size() <= 2) {
    return points;
  }

  PixelPoint lowest_point = FindLowestPoint(points);
  std::vector<PixelPoint> sorted_points = SortPointsByAngle(points, lowest_point);
  std::vector<PixelPoint> unique_points = RemoveCollinearPoints(sorted_points, lowest_point);

  if (unique_points.size() <= 1) {
    return HandleCollinearCase(unique_points, lowest_point);
  }

  return BuildHull(unique_points, lowest_point);
}

}  // namespace paramonov_v_bin_img_conv_hul_stl
