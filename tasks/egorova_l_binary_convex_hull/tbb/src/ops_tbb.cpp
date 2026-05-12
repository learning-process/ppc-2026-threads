#include "egorova_l_binary_convex_hull/tbb/include/ops_tbb.hpp"

#include <tbb/tbb.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <mutex>
#include <tuple>
#include <utility>
#include <vector>

#include "egorova_l_binary_convex_hull/common/include/common.hpp"
#include "oneapi/tbb/blocked_range2d.h"
#include "oneapi/tbb/parallel_for.h"

namespace egorova_l_binary_convex_hull {

namespace {

void ProcessNeighbours(const std::vector<uint8_t> &data, std::vector<int> &labels, std::vector<int> &queue, int cx,
                       int cy, int width, int height, int label) {
  const std::array<int, 4> dx = {1, -1, 0, 0};
  const std::array<int, 4> dy = {0, 0, 1, -1};
  for (int di = 0; di < 4; ++di) {
    const int nx = cx + dx.at(static_cast<std::size_t>(di));
    const int ny = cy + dy.at(static_cast<std::size_t>(di));
    if (nx < 0 || nx >= width || ny < 0 || ny >= height) {
      continue;
    }
    const int nidx = (ny * width) + nx;
    if (data[static_cast<std::size_t>(nidx)] != 0 && labels[nidx] == -1) {
      labels[nidx] = label;
      queue.push_back(nidx);
    }
  }
}

std::vector<int> LabelComponents(const std::vector<uint8_t> &data, int width, int height, int &num_components) {
  const int n = width * height;
  std::vector<int> labels(n, -1);
  num_components = 0;

  std::vector<int> queue;
  queue.reserve(static_cast<std::size_t>(n));

  for (int start = 0; start < n; ++start) {
    if (data[static_cast<std::size_t>(start)] == 0 || labels[start] != -1) {
      continue;
    }
    const int label = num_components++;
    labels[start] = label;
    queue.clear();
    queue.push_back(start);

    for (std::size_t qi = 0; qi < queue.size(); ++qi) {
      const int cur = queue[qi];
      const int cx = cur % width;
      const int cy = cur / width;
      ProcessNeighbours(data, labels, queue, cx, cy, width, height, label);
    }
  }
  return labels;
}

int64_t Cross(const Point &o, const Point &a, const Point &b) {
  return (static_cast<int64_t>(a.x - o.x) * (b.y - o.y)) - (static_cast<int64_t>(a.y - o.y) * (b.x - o.x));
}

std::vector<Point> ConvexHull(std::vector<Point> pts) {
  std::ranges::sort(pts,
                    [](const Point &lhs, const Point &rhs) { return std::tie(lhs.x, lhs.y) < std::tie(rhs.x, rhs.y); });
  pts.erase(std::ranges::unique(pts,
                                [](const Point &lhs, const Point &rhs) {
    return lhs.x == rhs.x && lhs.y == rhs.y;
  }).begin(),
            pts.end());

  if (pts.size() < 3) {
    return pts;
  }

  std::vector<Point> hull;
  hull.reserve(2 * pts.size());

  for (const auto &p : pts) {
    while (hull.size() >= 2 && Cross(hull[hull.size() - 2], hull.back(), p) <= 0) {
      hull.pop_back();
    }
    hull.push_back(p);
  }

  const std::size_t lower_size = hull.size() + 1;
  for (int i = static_cast<int>(pts.size()) - 2; i >= 0; --i) {
    while (hull.size() >= lower_size && Cross(hull[hull.size() - 2], hull.back(), pts[i]) <= 0) {
      hull.pop_back();
    }
    hull.push_back(pts[i]);
  }

  hull.pop_back();
  return hull;
}

void CollectPoints(const std::vector<int> &labels, int num_components, int width,
                   std::vector<std::vector<Point>> &component_points, std::vector<std::mutex> &mutexes,
                   const tbb::blocked_range2d<int> &r) {
  std::vector<std::vector<Point>> local(static_cast<std::size_t>(num_components));
  for (int row = r.rows().begin(); row < r.rows().end(); ++row) {
    for (int col = r.cols().begin(); col < r.cols().end(); ++col) {
      const int idx = (row * width) + col;
      const int lbl = labels[static_cast<std::size_t>(idx)];
      if (lbl >= 0) {
        local[static_cast<std::size_t>(lbl)].push_back({col, row});
      }
    }
  }
  for (int ci = 0; ci < num_components; ++ci) {
    auto &bucket = local[static_cast<std::size_t>(ci)];
    if (!bucket.empty()) {
      std::scoped_lock lock(mutexes[static_cast<std::size_t>(ci)]);
      auto &dst = component_points[static_cast<std::size_t>(ci)];
      dst.insert(dst.end(), bucket.begin(), bucket.end());
    }
  }
}

}  // namespace

BinaryConvexHullTBB::BinaryConvexHullTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool BinaryConvexHullTBB::ValidationImpl() {
  const auto &img = GetInput();
  if (img.width <= 0 || img.height <= 0) {
    return false;
  }
  const std::size_t expected = static_cast<std::size_t>(img.width) * static_cast<std::size_t>(img.height);
  return img.data.size() == expected;
}

bool BinaryConvexHullTBB::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}

bool BinaryConvexHullTBB::RunImpl() {
  const auto &img = GetInput();
  const int width = img.width;
  const int height = img.height;

  int num_components = 0;
  const std::vector<int> labels = LabelComponents(img.data, width, height, num_components);

  if (num_components == 0) {
    GetOutput().clear();
    return true;
  }

  std::vector<std::vector<Point>> component_points(static_cast<std::size_t>(num_components));
  std::vector<std::mutex> mutexes(static_cast<std::size_t>(num_components));

  tbb::parallel_for(tbb::blocked_range2d<int>(0, height, 0, width), [&](const tbb::blocked_range2d<int> &r) {
    CollectPoints(labels, num_components, width, component_points, mutexes, r);
  });

  std::vector<std::vector<Point>> hulls(static_cast<std::size_t>(num_components));

  tbb::parallel_for(0, num_components, [&](int ci) {
    auto &pts = component_points[static_cast<std::size_t>(ci)];
    if (!pts.empty()) {
      hulls[static_cast<std::size_t>(ci)] = ConvexHull(pts);
    }
  });

  GetOutput().clear();
  for (auto &h : hulls) {
    if (!h.empty()) {
      GetOutput().push_back(std::move(h));
    }
  }

  return true;
}

bool BinaryConvexHullTBB::PostProcessingImpl() {
  return !GetOutput().empty() || GetInput().data.empty() ||
         std::all_of(GetInput().data.begin(), GetInput().data.end(), [](uint8_t v) { return v == 0; });
}

}  // namespace egorova_l_binary_convex_hull
