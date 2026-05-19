#include "kamalagin_a_binary_image_convex_hull/all/include/ops_all.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <thread>
#include <utility>
#include <vector>

#include "kamalagin_a_binary_image_convex_hull/common/include/common.hpp"
#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/parallel_for.h"
#include "util/include/util.hpp"

namespace kamalagin_a_binary_image_convex_hull {

namespace {

constexpr std::array<std::pair<int, int>, 4> kFourNeighbors = {{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};

int64_t Cross(const Point &o, const Point &a, const Point &b) {
  return (static_cast<int64_t>(a.x - o.x) * (b.y - o.y)) - (static_cast<int64_t>(a.y - o.y) * (b.x - o.x));
}

int64_t DistSq(const Point &a, const Point &b) {
  const int dx = a.x - b.x;
  const int dy = a.y - b.y;
  return (static_cast<int64_t>(dx) * dx) + (static_cast<int64_t>(dy) * dy);
}

size_t GrahamFindPivot(std::vector<Point> &pts) {
  size_t pivot = 0;
  const size_t n = pts.size();
  for (size_t i = 1; i < n; ++i) {
    if (pts[i].y < pts[pivot].y || (pts[i].y == pts[pivot].y && pts[i].x < pts[pivot].x)) {
      pivot = i;
    }
  }
  return pivot;
}

void GrahamCollinearReduce(std::vector<Point> &pts) {
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

void GrahamScan(std::vector<Point> &pts, Hull &out) {
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

void GrahamHull(std::vector<Point> &pts, Hull &out) {
  out.clear();
  const size_t n = pts.size();
  if (n <= 1) {
    if (n == 1) {
      out.push_back(pts[0]);
    }
    return;
  }
  const size_t pivot = GrahamFindPivot(pts);
  std::swap(pts[0], pts[pivot]);
  const Point &p0 = pts[0];
  std::sort(pts.begin() + 1, pts.end(), [&p0](const Point &a, const Point &b) {
    const int64_t c = Cross(p0, a, b);
    if (c != 0) {
      return c > 0;
    }
    return DistSq(p0, a) < DistSq(p0, b);
  });
  GrahamCollinearReduce(pts);
  GrahamScan(pts, out);
}

void FloodFillComponent(const BinaryImage &img, int start_row, int start_col, std::vector<int> &label,
                        std::vector<Point> &component_pts) {
  const int rows = img.rows;
  const int cols = img.cols;
  component_pts.clear();
  std::vector<std::pair<int, int>> stack;
  stack.emplace_back(start_row, start_col);
  const size_t start_idx = img.Index(start_row, start_col);
  label[start_idx] = 1;
  while (!stack.empty()) {
    const auto [cur_r, cur_c] = stack.back();
    stack.pop_back();
    component_pts.push_back(Point{.x = cur_c, .y = cur_r});
    for (const auto &[dr, dc] : kFourNeighbors) {
      const int nr = cur_r + dr;
      const int nc = cur_c + dc;
      if (nr < 0 || nr >= rows || nc < 0 || nc >= cols) {
        continue;
      }
      const size_t nidx = (static_cast<size_t>(nr) * static_cast<size_t>(cols)) + static_cast<size_t>(nc);
      if (img.data[nidx] != 0 && label[nidx] == 0) {
        label[nidx] = 1;
        stack.emplace_back(nr, nc);
      }
    }
  }
}

void CollectComponents(const BinaryImage &img, std::vector<std::vector<Point>> &components) {
  components.clear();
  const int rows = img.rows;
  const int cols = img.cols;
  const size_t total = static_cast<size_t>(rows) * static_cast<size_t>(cols);
  std::vector<int> label(total, 0);
  for (int row = 0; row < rows; ++row) {
    for (int col = 0; col < cols; ++col) {
      const size_t idx = img.Index(row, col);
      if (img.data[idx] == 0 || label[idx] != 0) {
        continue;
      }
      components.emplace_back();
      FloodFillComponent(img, row, col, label, components.back());
    }
  }
}

void RunHullRangeOmp(std::vector<std::vector<Point>> &components, HullList &hulls, size_t range_begin,
                     size_t range_end) {
  const int n = static_cast<int>(range_end - range_begin);
  if (n <= 0) {
    return;
  }
#pragma omp parallel for default(none) shared(components, hulls, range_begin, n) schedule(dynamic)
  for (int i = 0; i < n; ++i) {
    const auto idx = range_begin + static_cast<size_t>(i);
    GrahamHull(components[idx], hulls[idx]);
  }
}

void RunHullRangeStl(std::vector<std::vector<Point>> &components, HullList &hulls, size_t range_begin,
                     size_t range_end) {
  const size_t len = range_end - range_begin;
  if (len == 0) {
    return;
  }
  const auto requested = static_cast<size_t>(ppc::util::GetNumThreads());
  const size_t num_threads = std::min<size_t>(std::max<size_t>(requested, 1), len);
  const size_t chunk = len / num_threads;
  const size_t remainder = len % num_threads;
  std::vector<std::thread> threads;
  threads.reserve(num_threads);
  size_t offset = 0;
  for (size_t ti = 0; ti < num_threads; ++ti) {
    const size_t span = chunk + (ti < remainder ? 1 : 0);
    const size_t local_begin = range_begin + offset;
    const size_t local_end = local_begin + span;
    offset += span;
    threads.emplace_back([&components, &hulls, local_begin, local_end]() {
      for (size_t i = local_begin; i < local_end; ++i) {
        GrahamHull(components[i], hulls[i]);
      }
    });
  }
  for (auto &th : threads) {
    th.join();
  }
}

void RunHullRangeTbb(std::vector<std::vector<Point>> &components, HullList &hulls, size_t range_begin,
                     size_t range_end) {
  if (range_begin >= range_end) {
    return;
  }
  tbb::parallel_for(tbb::blocked_range<size_t>(range_begin, range_end), [&](const tbb::blocked_range<size_t> &range) {
    for (size_t i = range.begin(); i != range.end(); ++i) {
      GrahamHull(components[i], hulls[i]);
    }
  });
}

void RunBinaryImageConvexHullAll(const BinaryImage &img, HullList &hulls) {
  hulls.clear();
  std::vector<std::vector<Point>> components;
  CollectComponents(img, components);
  const size_t count = components.size();
  hulls.resize(count);
  if (count == 0) {
    return;
  }
  const size_t b1 = count / 3;
  const size_t b2 = (2 * count) / 3;
  RunHullRangeOmp(components, hulls, 0, b1);
  RunHullRangeStl(components, hulls, b1, b2);
  RunHullRangeTbb(components, hulls, b2, count);
}

}  // namespace

KamalaginABinaryImageConvexHullALL::KamalaginABinaryImageConvexHullALL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = HullList{};
}

bool KamalaginABinaryImageConvexHullALL::ValidationImpl() {
  const auto &img = GetInput();
  if (img.rows < 0 || img.cols < 0) {
    return false;
  }
  if (img.rows == 0 || img.cols == 0) {
    return img.data.empty();
  }
  if (img.rows > 8192 || img.cols > 8192) {
    return false;
  }
  return (static_cast<size_t>(img.rows) * static_cast<size_t>(img.cols)) == img.data.size();
}

bool KamalaginABinaryImageConvexHullALL::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}

bool KamalaginABinaryImageConvexHullALL::RunImpl() {
  RunBinaryImageConvexHullAll(GetInput(), GetOutput());
  return true;
}

bool KamalaginABinaryImageConvexHullALL::PostProcessingImpl() {
  return true;
}

}  // namespace kamalagin_a_binary_image_convex_hull
