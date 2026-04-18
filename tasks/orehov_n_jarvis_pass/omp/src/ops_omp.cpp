#include "orehov_n_jarvis_pass/omp/include/ops_omp.hpp"

#include <omp.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include "orehov_n_jarvis_pass/common/include/common.hpp"

namespace orehov_n_jarvis_pass {

OrehovNJarvisPassOMP::OrehovNJarvisPassOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<Point>();
}

bool OrehovNJarvisPassOMP::ValidationImpl() {
  return !GetInput().empty();
}

bool OrehovNJarvisPassOMP::PreProcessingImpl() {
  return true;
}

bool OrehovNJarvisPassOMP::RunImpl() {
  const auto &input = GetInput();

  if (input.size() == 1 || input.size() == 2) {
    GetOutput() = input;
    return true;
  }

  Point current = FindFirstElem(input);
  GetOutput().push_back(current);

  while (true) {
    Point next = FindNext(current, input);
    if (next == GetOutput()[0]) {
      break;
    }

    current = next;
    GetOutput().push_back(next);
  }

  return true;
}

Point OrehovNJarvisPassOMP::FindNext(Point current, const std::vector<Point> &input) const {
  const size_t n = input.size();
  Point initial_candidate = (current == input[0]) ? input[1] : input[0];

  int max_threads = omp_get_max_threads();
  std::vector<Point> local_bests(max_threads, initial_candidate);
  const int n_int = static_cast<int>(n);

#pragma omp parallel default(none) shared(input, n_int, current, local_bests)
  {
    int tid = omp_get_thread_num();
    Point thread_local_best = local_bests[tid];

#pragma omp for nowait
    for (int i = 0; i < n_int; ++i) {
      const Point &point = input[i];
      if (current == point) {
        continue;
      }

      double orient = CheckLeft(current, thread_local_best, point);

      if (orient > 0) {
        thread_local_best = point;
      } else if (std::abs(orient) < 1e-9) {
        if (DistanceSquared(current, point) > DistanceSquared(current, thread_local_best)) {
          thread_local_best = point;
        }
      }
    }
    local_bests[tid] = thread_local_best;
  }

  Point global_next = local_bests[0];
  for (int i = 1; i < max_threads; ++i) {
    double orient = CheckLeft(current, global_next, local_bests[i]);
    if (orient > 0) {
      global_next = local_bests[i];
    } else if (std::abs(orient) < 1e-9) {
      if (DistanceSquared(current, local_bests[i]) > DistanceSquared(current, global_next)) {
        global_next = local_bests[i];
      }
    }
  }

  return global_next;
}

double OrehovNJarvisPassOMP::CheckLeft(Point a, Point b, Point c) {
  return ((b.x - a.x) * (c.y - a.y)) - ((b.y - a.y) * (c.x - a.x));
}

Point OrehovNJarvisPassOMP::FindFirstElem(const std::vector<Point> &input) const {
  Point current = input[0];
  for (const auto &f : input) {
    if (f.x < current.x || (f.y < current.y && f.x == current.x)) {
      current = f;
    }
  }
  return current;
}

double OrehovNJarvisPassOMP::DistanceSquared(Point a, Point b) {
  double dx = a.x - b.x;
  double dy = a.y - b.y;
  return dx * dx + dy * dy;
}

bool OrehovNJarvisPassOMP::PostProcessingImpl() {
  return true;
}

}  // namespace orehov_n_jarvis_pass
