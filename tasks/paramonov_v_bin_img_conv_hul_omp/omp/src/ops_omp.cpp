#include "paramonov_v_bin_img_conv_hul_omp/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <stack>
#include <utility>
#include <vector>

#include "paramonov_v_bin_img_conv_hul_omp/common/include/common.hpp"

namespace paramonov_v_bin_img_conv_hul_omp {

namespace {
constexpr std::array<std::pair<int, int>, 4> kNeighbors = {{{1, 0}, {-1, 0}, {0, 1}, {0, -1}}};

bool ComparePoints(const PixelPoint &a, const PixelPoint &b) {
  if (a.row != b.row) {
    return a.row < b.row;
  }
  return a.col < b.col;
}

}  // namespace

ConvexHullOMP::ConvexHullOMP(const InputType &input) {
  SetTypeOfTask(StaticTaskType());
  GetInput() = input;
}

bool ConvexHullOMP::ValidationImpl() {
  const auto &img = GetInput();
  if (img.rows <= 0 || img.cols <= 0) {
    return false;
  }

  const size_t expected_size = static_cast<size_t>(img.rows) * static_cast<size_t>(img.cols);
  return img.pixels.size() == expected_size;
}

bool ConvexHullOMP::PreProcessingImpl() {
  working_image_ = GetInput();
  BinarizeImage();
  GetOutput().clear();
  return true;
}

bool ConvexHullOMP::RunImpl() {
  ExtractConnectedComponents();
  return true;
}

bool ConvexHullOMP::PostProcessingImpl() {
  return true;
}

void ConvexHullOMP::BinarizeImage(uint8_t threshold) {
  const size_t size = working_image_.pixels.size();
  auto &pixels = working_image_.pixels;

#pragma omp parallel for default(none) shared(pixels, threshold, size)
  for (size_t i = 0; i < size; ++i) {
    pixels[i] = pixels[i] > threshold ? uint8_t{255} : uint8_t{0};
  }
}

void ConvexHullOMP::FloodFill(int start_row, int start_col, std::vector<bool> &visited,
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

    for (const auto &neighbor : kNeighbors) {
      int dr = neighbor.first;
      int dc = neighbor.second;
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

void ConvexHullOMP::FindStartPoints(const std::vector<uint8_t> &pixels, int rows, int cols, std::vector<bool> &visited,
                                    std::vector<std::pair<int, int>> &start_points) {
#pragma omp parallel for default(none) shared(rows, cols, pixels, visited, start_points)
  for (int row = 0; row < rows; ++row) {
    for (int col = 0; col < cols; ++col) {
      size_t idx = (static_cast<size_t>(row) * static_cast<size_t>(cols)) + static_cast<size_t>(col);
      if (pixels[idx] == 255 && !visited[idx]) {
#pragma omp critical
        {
          if (!visited[idx]) {
            visited[idx] = true;
            start_points.emplace_back(row, col);
          }
        }
      }
    }
  }
}

void ConvexHullOMP::ProcessComponent(int start_row, int start_col, int rows, int cols, size_t total_pixels,
                                     const std::vector<uint8_t> &pixels,
                                     std::vector<std::vector<PixelPoint>> &components) {
  std::vector<bool> local_visited(total_pixels, false);
  std::vector<PixelPoint> component;

  std::stack<PixelPoint> pixel_stack;
  pixel_stack.emplace(start_row, start_col);
  local_visited[(static_cast<size_t>(start_row) * static_cast<size_t>(cols)) + static_cast<size_t>(start_col)] = true;

  while (!pixel_stack.empty()) {
    PixelPoint current = pixel_stack.top();
    pixel_stack.pop();
    component.push_back(current);

    for (const auto &neighbor : kNeighbors) {
      int dr = neighbor.first;
      int dc = neighbor.second;
      int next_row = current.row + dr;
      int next_col = current.col + dc;

      if (next_row >= 0 && next_row < rows && next_col >= 0 && next_col < cols) {
        size_t idx = (static_cast<size_t>(next_row) * static_cast<size_t>(cols)) + static_cast<size_t>(next_col);
        if (!local_visited[idx] && pixels[idx] == 255) {
          local_visited[idx] = true;
          pixel_stack.emplace(next_row, next_col);
        }
      }
    }
  }

  if (!component.empty()) {
    std::vector<PixelPoint> hull = ComputeConvexHull(component);
#pragma omp critical
    {
      components.push_back(std::move(hull));
    }
  }
}

void ConvexHullOMP::ExtractConnectedComponents() {
  const int rows = working_image_.rows;
  const int cols = working_image_.cols;
  const size_t total_pixels = static_cast<size_t>(rows) * static_cast<size_t>(cols);

  std::vector<bool> visited(total_pixels, false);
  std::vector<std::vector<PixelPoint>> components;
  std::vector<std::pair<int, int>> start_points;

  auto &pixels = working_image_.pixels;

  FindStartPoints(pixels, rows, cols, visited, start_points);

  const size_t num_points = start_points.size();
#pragma omp parallel for default(none) shared(num_points, start_points, total_pixels, rows, cols, pixels, components)
  for (size_t i = 0; i < num_points; ++i) {
    const auto &start_point = start_points[i];
    int start_row = start_point.first;
    int start_col = start_point.second;
    ProcessComponent(start_row, start_col, rows, cols, total_pixels, pixels, components);
  }

  GetOutput() = std::move(components);
}

int64_t ConvexHullOMP::Orientation(const PixelPoint &p, const PixelPoint &q, const PixelPoint &r) {
  return (static_cast<int64_t>(q.col - p.col) * (r.row - p.row)) -
         (static_cast<int64_t>(q.row - p.row) * (r.col - p.col));
}

std::vector<PixelPoint> ConvexHullOMP::ComputeConvexHull(const std::vector<PixelPoint> &points) {
  if (points.size() <= 2) {
    return points;
  }

  auto lowest_point = *std::ranges::min_element(points, ComparePoints);

  std::vector<PixelPoint> sorted_points;
  sorted_points.reserve(points.size() - 1);
  std::ranges::copy_if(points, std::back_inserter(sorted_points), [&lowest_point](const PixelPoint &p) {
    return (p.row != lowest_point.row) || (p.col != lowest_point.col);
  });

  std::ranges::sort(sorted_points, [&lowest_point](const PixelPoint &a, const PixelPoint &b) {
    int64_t orient = Orientation(lowest_point, a, b);
    if (orient == 0) {
      int64_t dist_a = ((a.row - lowest_point.row) * (a.row - lowest_point.row)) +
                       ((a.col - lowest_point.col) * (a.col - lowest_point.col));
      int64_t dist_b = ((b.row - lowest_point.row) * (b.row - lowest_point.row)) +
                       ((b.col - lowest_point.col) * (b.col - lowest_point.col));
      return dist_a < dist_b;
    }
    return orient > 0;
  });

  std::vector<PixelPoint> hull;
  hull.reserve(points.size());
  hull.push_back(lowest_point);

  for (const auto &p : sorted_points) {
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

}  // namespace paramonov_v_bin_img_conv_hul_omp
