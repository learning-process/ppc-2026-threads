#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "kamalagin_a_binary_image_convex_hull/common/include/common.hpp"

namespace kamalagin_a_binary_image_convex_hull::detail {

constexpr std::array<std::pair<int, int>, 4> kFourNeighbors = {{
    {-1, 0},
    {1, 0},
    {0, -1},
    {0, 1},
}};

inline int64_t Cross(const Point &o, const Point &a, const Point &b) {
  return (static_cast<int64_t>(a.x - o.x) * (b.y - o.y)) - (static_cast<int64_t>(a.y - o.y) * (b.x - o.x));
}

inline int64_t DistSq(const Point &a, const Point &b) {
  int dx = a.x - b.x;
  int dy = a.y - b.y;
  return (static_cast<int64_t>(dx) * dx) + (static_cast<int64_t>(dy) * dy);
}

inline size_t GrahamFindPivot(std::vector<Point> &pts) {
  size_t pivot = 0;
  const size_t n = pts.size();
  for (size_t i = 1; i < n; ++i) {
    if (pts[i].y < pts[pivot].y || (pts[i].y == pts[pivot].y && pts[i].x < pts[pivot].x)) {
      pivot = i;
    }
  }
  return pivot;
}

inline void GrahamCollinearReduce(std::vector<Point> &pts) {
  size_t m = 1;
  const size_t n = pts.size();
  for (size_t i = 2; i < n; ++i) {
    while (m > 0 && Cross(pts[m - 1], pts[m], pts[i]) == 0) {
      --m;
    }
    ++m;
    pts[m] = pts[i];
  }
  pts.resize(m + 1);
}

inline void GrahamScan(std::vector<Point> &pts, Hull &out) {
  out.push_back(pts[0]);
  if (pts.size() <= 2) {
    if (pts.size() == 2) {
      out.push_back(pts[1]);
    }
    return;
  }
  out.push_back(pts[1]);
  for (size_t i = 2; i < pts.size(); ++i) {
    while (out.size() >= 2 && Cross(out[out.size() - 2], out.back(), pts[i]) <= 0) {
      out.pop_back();
    }
    out.push_back(pts[i]);
  }
}

inline void GrahamHull(std::vector<Point> &pts, Hull &out) {
  out.clear();
  const size_t n = pts.size();
  if (n <= 1) {
    if (n == 1) {
      out.push_back(pts[0]);
    }
    return;
  }
  size_t pivot = GrahamFindPivot(pts);
  std::swap(pts[0], pts[pivot]);
  const Point &p0 = pts[0];
  std::sort(pts.begin() + 1, pts.end(), [&p0](const Point &a, const Point &b) {
    int64_t c = Cross(p0, a, b);
    if (c != 0) {
      return c > 0;
    }
    return DistSq(p0, a) < DistSq(p0, b);
  });
  GrahamCollinearReduce(pts);
  GrahamScan(pts, out);
}

inline void FloodFillComponent(const BinaryImage &img, int start_row, int start_col, std::vector<int> &label,
                               std::vector<Point> &component_pts) {
  const int rows = img.rows;
  const int cols = img.cols;
  component_pts.clear();
  std::vector<std::pair<int, int>> stack;
  stack.emplace_back(start_row, start_col);
  const size_t start_idx = img.Index(start_row, start_col);
  label[start_idx] = 1;
  while (!stack.empty()) {
    auto [cur_r, cur_c] = stack.back();
    stack.pop_back();
    component_pts.push_back(Point{.x = cur_c, .y = cur_r});
    for (const auto &[dr, dc] : kFourNeighbors) {
      int nr = cur_r + dr;
      int nc = cur_c + dc;
      if (nr < 0 || nr >= rows || nc < 0 || nc >= cols) {
        continue;
      }
      size_t nidx = (static_cast<size_t>(nr) * static_cast<size_t>(cols)) + static_cast<size_t>(nc);
      if (img.data[nidx] != 0 && label[nidx] == 0) {
        label[nidx] = 1;
        stack.emplace_back(nr, nc);
      }
    }
  }
}

/// Sequential: 4-connected components + Graham hull per component.
inline void RunBinaryImageConvexHull(const BinaryImage &img, HullList &hulls) {
  hulls.clear();
  const int rows = img.rows;
  const int cols = img.cols;
  const size_t total = static_cast<size_t>(rows) * static_cast<size_t>(cols);
  std::vector<int> label(total, 0);
  std::vector<Point> component_pts;
  component_pts.reserve(total);
  for (int row = 0; row < rows; ++row) {
    for (int col = 0; col < cols; ++col) {
      size_t idx = img.Index(row, col);
      if (img.data[idx] == 0 || label[idx] != 0) {
        continue;
      }
      FloodFillComponent(img, row, col, label, component_pts);
      Hull hull;
      GrahamHull(component_pts, hull);
      hulls.push_back(std::move(hull));
    }
  }
}

}  // namespace kamalagin_a_binary_image_convex_hull::detail
