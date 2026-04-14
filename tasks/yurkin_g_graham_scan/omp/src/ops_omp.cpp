#include "yurkin_g_graham_scan/omp/include/ops_omp.hpp"

#include <algorithm>
#include <cstddef>
#include <ranges>
#include <vector>

#if defined(__cpp_lib_execution) && (__cpp_lib_execution >= 201603)
#  include <execution>
#endif

#include "yurkin_g_graham_scan/common/include/common.hpp"

namespace yurkin_g_graham_scan {

YurkinGGrahamScanOMP::YurkinGGrahamScanOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput().clear();
}

bool YurkinGGrahamScanOMP::ValidationImpl() {
  return !GetInput().empty();
}

bool YurkinGGrahamScanOMP::PreProcessingImpl() {
  auto &pts = GetInput();
  if (pts.empty()) {
    return true;
  }

  // Сортировка по x, затем по y. При наличии поддержки execution policy
  // используем параллельную политику для ускорения.
#if defined(__cpp_lib_execution) && (__cpp_lib_execution >= 201603)
  std::sort(std::execution::par, pts.begin(), pts.end(), [](const Point &a, const Point &b) {
    if (a.x != b.x) {
      return a.x < b.x;
    }
    return a.y < b.y;
  });
#else
  std::sort(pts.begin(), pts.end(), [](const Point &a, const Point &b) {
    if (a.x != b.x) {
      return a.x < b.x;
    }
    return a.y < b.y;
  });
#endif

  // Удаление дубликатов (последовательный проход)
  std::vector<Point> tmp;
  tmp.reserve(pts.size());
  for (const auto &p : pts) {
    if (tmp.empty() || tmp.back().x != p.x || tmp.back().y != p.y) {
      tmp.push_back(p);
    }
  }
  pts.swap(tmp);

  return !pts.empty();
}

namespace {  // helper: cross product
long double Cross(const Point &o, const Point &a, const Point &b) {
  return (static_cast<long double>(a.x - o.x) * static_cast<long double>(b.y - o.y)) -
         (static_cast<long double>(a.y - o.y) * static_cast<long double>(b.x - o.x));
}
}  // namespace

bool YurkinGGrahamScanOMP::RunImpl() {
  const InType pts_in = GetInput();
  const std::size_t n = pts_in.size();
  if (n == 0) {
    GetOutput().clear();
    return true;
  }
  if (n == 1) {
    GetOutput() = pts_in;
    return true;
  }

  // Локальная копия точек; сортировка для надёжности (PreProcessing уже сортировал)
  InType pts = pts_in;
#if defined(__cpp_lib_execution) && (__cpp_lib_execution >= 201603)
  std::sort(std::execution::par, pts.begin(), pts.end(), [](const Point &a, const Point &b) {
    if (a.x != b.x) {
      return a.x < b.x;
    }
    return a.y < b.y;
  });
#else
  std::sort(pts.begin(), pts.end(), [](const Point &a, const Point &b) {
    if (a.x != b.x) {
      return a.x < b.x;
    }
    return a.y < b.y;
  });
#endif

  // Построение нижней оболочки
  OutType lower;
  lower.reserve(pts.size());
  for (const auto &p : pts) {
    while (lower.size() >= 2 && Cross(lower[lower.size() - 2], lower[lower.size() - 1], p) <= 0) {
      lower.pop_back();
    }
    lower.push_back(p);
  }

  // Построение верхней оболочки (range-based обратный проход)
  OutType upper;
  upper.reserve(pts.size());
  for (const auto &p : std::ranges::reverse_view(pts)) {
    while (upper.size() >= 2 && Cross(upper[upper.size() - 2], upper[upper.size() - 1], p) <= 0) {
      upper.pop_back();
    }
    upper.push_back(p);
  }

  // Объединение (без дублирования крайних точек)
  OutType hull;
  hull.reserve(lower.size() + upper.size());
  for (const auto &pt : lower) {
    hull.push_back(pt);
  }
  for (std::size_t i = 1; i + 1 < upper.size(); ++i) {
    hull.push_back(upper[i]);
  }

  GetOutput() = hull;
  return true;
}

bool YurkinGGrahamScanOMP::PostProcessingImpl() {
  if (GetInput().empty()) {
    return true;
  }
  return !GetOutput().empty();
}

}  // namespace yurkin_g_graham_scan
