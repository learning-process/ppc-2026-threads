#include "orehov_n_jarvis_pass/stl/include/ops_stl.hpp"

#include <cmath>
#include <cstddef>
#include <set>
#include <thread>
#include <vector>

#include "orehov_n_jarvis_pass/common/include/common.hpp"

namespace orehov_n_jarvis_pass {

namespace {

struct BestState {
  Point point;
  bool valid = false;
};

bool IsBetterPoint(Point current, Point candidate, Point best) {
  double orient = OrehovNJarvisPassSTL::CheckLeft(current, best, candidate);
  if (orient > 0.0) {
    return true;
  }
  if (orient == 0.0) {
    double dx_c = candidate.x - current.x;
    double dy_c = candidate.y - current.y;
    double dist_c = (dx_c * dx_c) + (dy_c * dy_c);

    double dx_b = best.x - current.x;
    double dy_b = best.y - current.y;
    double dist_b = (dx_b * dx_b) + (dy_b * dy_b);

    return dist_c > dist_b;
  }
  return false;
}

BestState FindBestInRange(const std::vector<Point> &points, Point current, size_t start, size_t end) {
  BestState local{.point = Point(), .valid = false};
  for (size_t j = start; j < end; ++j) {
    const Point &p = points[j];
    if (p == current) {
      continue;
    }
    if (!local.valid || IsBetterPoint(current, p, local.point)) {
      local.point = p;
      local.valid = true;
    }
  }
  return local;
}

BestState CombineBestStates(const BestState &a, const BestState &b, Point current) {
  if (!a.valid) {
    return b;
  }
  if (!b.valid) {
    return a;
  }
  return BestState{.point = IsBetterPoint(current, b.point, a.point) ? b.point : a.point, .valid = true};
}

}  // namespace

OrehovNJarvisPassSTL::OrehovNJarvisPassSTL(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<Point>();
}

bool OrehovNJarvisPassSTL::ValidationImpl() {
  return (!GetInput().empty());
}

bool OrehovNJarvisPassSTL::PreProcessingImpl() {
  std::set<Point> tmp(GetInput().begin(), GetInput().end());
  input_.assign(tmp.begin(), tmp.end());
  return true;
}

bool OrehovNJarvisPassSTL::RunImpl() {
  if (input_.size() == 1 || input_.size() == 2) {
    res_ = input_;
    return true;
  }

  Point current = FindFirstElem();
  res_.push_back(current);

  while (true) {
    Point next = FindNext(current);
    if (next == res_[0]) {
      break;
    }
    current = next;
    res_.push_back(next);
  }

  return true;
}

Point OrehovNJarvisPassSTL::FindNext(Point current) const {
  const size_t n = input_.size();
  if (n == 0) {
    return current;
  }

  unsigned int num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 1;
  }
  if (num_threads > n) {
    num_threads = static_cast<unsigned int>(n);
  }

  size_t chunk = n / num_threads;
  size_t remainder = n % num_threads;
  size_t start = 0;

  std::vector<std::thread> threads;
  std::vector<BestState> results(num_threads, BestState{.point = Point(), .valid = false});

  for (unsigned int i = 0; i < num_threads; ++i) {
    size_t end = start + chunk + (i < remainder ? 1 : 0);
    threads.emplace_back(
        [this, &results, i, start, end, current]() { results[i] = FindBestInRange(input_, current, start, end); });
    start = end;
  }

  for (auto &t : threads) {
    t.join();
  }

  BestState final_best{.point = Point(), .valid = false};
  for (const auto &best : results) {
    final_best = CombineBestStates(final_best, best, current);
  }

  return final_best.valid ? final_best.point : current;
}

double OrehovNJarvisPassSTL::CheckLeft(Point a, Point b, Point c) {
  return ((b.x - a.x) * (c.y - a.y)) - ((b.y - a.y) * (c.x - a.x));
}

Point OrehovNJarvisPassSTL::FindFirstElem() const {
  Point current = input_[0];
  for (const auto &f : input_) {
    if (f.x < current.x || (f.y < current.y && f.x == current.x)) {
      current = f;
    }
  }
  return current;
}

double OrehovNJarvisPassSTL::Distance(Point a, Point b) {
  return std::sqrt(std::pow(a.y - b.y, 2) + std::pow(a.x - b.x, 2));
}

bool OrehovNJarvisPassSTL::PostProcessingImpl() {
  GetOutput() = res_;
  return true;
}

}  // namespace orehov_n_jarvis_pass
