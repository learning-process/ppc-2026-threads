#include "egorova_l_binary_convex_hull/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <queue>
#include <utility>
#include <vector>

#include "egorova_l_binary_convex_hull/common/include/common.hpp"

namespace egorova_l_binary_convex_hull {

namespace {
struct Pixel {
  int x;
  int y;
  Pixel(int x_coord, int y_coord) : x(x_coord), y(y_coord) {}
};

// Вспомогательная функция для проверки соседей
inline void CheckNeighbor(const std::vector<uint8_t> &image, int width, int height, int next_x, int next_y, int label,
                          std::vector<int> &labels, std::queue<Pixel> &queue) {
  if (next_x >= 0 && next_x < width && next_y >= 0 && next_y < height) {
    const size_t next_index = (static_cast<size_t>(next_y) * static_cast<size_t>(width)) + static_cast<size_t>(next_x);
    if (image[next_index] != 0 && labels[next_index] == 0) {
      labels[next_index] = label;
      queue.emplace(next_x, next_y);
    }
  }
}
}  // namespace

// ==================== Публичные методы ====================

BinaryConvexHullOMP::BinaryConvexHullOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
}

bool BinaryConvexHullOMP::ValidationImpl() {
  const auto &input = GetInput();
  return input.width > 0 && input.height > 0 && !input.data.empty() &&
         input.data.size() == static_cast<size_t>(input.width) * static_cast<size_t>(input.height);
}

bool BinaryConvexHullOMP::PreProcessingImpl() {
  GetOutput().clear();
  return true;
}

bool BinaryConvexHullOMP::RunImpl() {
  const int width = GetInput().width;
  const int height = GetInput().height;
  const auto &image = GetInput().data;

  std::vector<std::vector<Point>> components = FindComponents(image, width, height);

  auto &output = GetOutput();
  output.resize(components.size());

  omp_set_num_threads(omp_get_max_threads());

  const int num_components = static_cast<int>(components.size());
#pragma omp parallel for schedule(dynamic, 1) default(none) shared(components, output, num_components)
  for (int i = 0; i < num_components; ++i) {
    if (components[i].empty()) {
      continue;
    }

    std::vector<Point> hull;
    BuildConvexHull(components[i], hull);

#pragma omp critical
    {
      output[i] = std::move(hull);
    }
  }

  return true;
}

bool BinaryConvexHullOMP::PostProcessingImpl() {
  return true;
}

// ==================== Приватные методы ====================

std::vector<std::vector<Point>> BinaryConvexHullOMP::FindComponents(const std::vector<uint8_t> &image, int width,
                                                                    int height) {
  const size_t image_size = static_cast<size_t>(width) * static_cast<size_t>(height);
  std::vector<int> labels(image_size, 0);
  int label_counter = 0;

  std::vector<std::vector<Point>> components;
  components.reserve(100);

  for (int row = 0; row < height; ++row) {
    for (int col = 0; col < width; ++col) {
      const size_t index = (static_cast<size_t>(row) * static_cast<size_t>(width)) + static_cast<size_t>(col);

      if (image[index] != 0 && labels[index] == 0) {
        ++label_counter;
        std::vector<Point> component_points;

        ProcessComponent(image, width, height, col, row, label_counter, labels, component_points);

        if (!component_points.empty()) {
          components.push_back(std::move(component_points));
        }
      }
    }
  }

  return components;
}

void BinaryConvexHullOMP::ProcessComponent(const std::vector<uint8_t> &image, int width, int height, int start_x,
                                           int start_y, int label, std::vector<int> &labels,
                                           std::vector<Point> &component_points) {
  if (start_x < 0 || start_x >= width || start_y < 0 || start_y >= height) {
    component_points.clear();
    return;
  }

  std::queue<Pixel> queue;
  queue.emplace(start_x, start_y);

  const size_t start_index = (static_cast<size_t>(start_y) * static_cast<size_t>(width)) + static_cast<size_t>(start_x);
  labels[start_index] = label;

  component_points.clear();
  component_points.reserve(1000);

  while (!queue.empty()) {
    Pixel current = queue.front();
    queue.pop();

    component_points.emplace_back(current.x, current.y);

    // Явно перебираем все направления
    CheckNeighbor(image, width, height, current.x + 1, current.y, label, labels, queue);
    CheckNeighbor(image, width, height, current.x - 1, current.y, label, labels, queue);
    CheckNeighbor(image, width, height, current.x, current.y + 1, label, labels, queue);
    CheckNeighbor(image, width, height, current.x, current.y - 1, label, labels, queue);
  }
}

int64_t BinaryConvexHullOMP::CrossProduct(const Point &a, const Point &b, const Point &c) {
  const auto dx1 = static_cast<int64_t>(b.x - a.x);
  const auto dy1 = static_cast<int64_t>(c.y - a.y);
  const auto dx2 = static_cast<int64_t>(b.y - a.y);
  const auto dy2 = static_cast<int64_t>(c.x - a.x);
  return (dx1 * dy1) - (dx2 * dy2);
}

void BinaryConvexHullOMP::BuildConvexHull(std::vector<Point> &points, std::vector<Point> &hull) {
  const size_t points_count = points.size();

  if (points_count == 0) {
    hull.clear();
    return;
  }

  if (points_count <= 2) {
    hull = points;
    return;
  }

  std::ranges::sort(
      points, [](const Point &lhs, const Point &rhs) { return lhs.x < rhs.x || (lhs.x == rhs.x && lhs.y < rhs.y); });

  // Нижняя оболочка
  hull.clear();
  hull.reserve(points_count + 1);
  for (size_t i = 0; i < points_count; ++i) {
    while (hull.size() >= 2 && CrossProduct(hull[hull.size() - 2], hull.back(), points[i]) <= 0) {
      hull.pop_back();
    }
    hull.push_back(points[i]);
  }

  // Верхняя оболочка
  const size_t lower_hull_size = hull.size();
  for (size_t i = points_count - 1; i > 0; --i) {
    const size_t idx = i - 1;
    while (hull.size() > lower_hull_size && CrossProduct(hull[hull.size() - 2], hull.back(), points[idx]) <= 0) {
      hull.pop_back();
    }
    hull.push_back(points[idx]);
  }

  if (hull.size() > 1 && hull.front().x == hull.back().x && hull.front().y == hull.back().y) {
    hull.pop_back();
  }
}

}  // namespace egorova_l_binary_convex_hull
