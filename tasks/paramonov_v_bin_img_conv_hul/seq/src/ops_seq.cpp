#include "paramonov_v_bin_img_conv_hul/seq/include/ops_seq.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iterator>  // для std::back_inserter
#include <stack>
#include <utility>  // для std::pair
#include <vector>

namespace paramonov_v_bin_img_conv_hul {

namespace {
constexpr std::array<std::pair<int, int>, 4> kNeighbors = {{{1, 0}, {-1, 0}, {0, 1}, {0, -1}}};

bool ComparePoints(const PixelPoint &a, const PixelPoint &b) {
  if (a.row != b.row) {
    return a.row < b.row;
  }
  return a.col < b.col;
}

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
  std::ranges::transform(working_image_.pixels, working_image_.pixels.begin(),
                         [threshold](uint8_t pixel) { return pixel > threshold ? uint8_t{255} : uint8_t{0}; });
}

void ConvexHullSequential::FloodFill(int start_row, int start_col, std::vector<bool> &visited,
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

void ConvexHullSequential::ExtractConnectedComponents() {
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

int64_t ConvexHullSequential::Orientation(const PixelPoint &p, const PixelPoint &q, const PixelPoint &r) {
  return (static_cast<int64_t>(q.col - p.col) * (r.row - p.row)) -
         (static_cast<int64_t>(q.row - p.row) * (r.col - p.col));
}

std::vector<PixelPoint> ConvexHullSequential::ComputeConvexHull(const std::vector<PixelPoint> &points) const {
  if (points.size() <= 2) {
    return points;
  }

  // Находим точку с наименьшими координатами
  auto lowest_point = *std::ranges::min_element(points, ComparePoints);

  // Копируем и сортируем по полярному углу
  std::vector<PixelPoint> sorted_points;
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

  // Строим выпуклую оболочку
  std::vector<PixelPoint> hull;
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

}  // namespace paramonov_v_bin_img_conv_hul
