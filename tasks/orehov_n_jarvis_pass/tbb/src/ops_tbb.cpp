#include "orehov_n_jarvis_pass/tbb/include/ops_tbb.hpp"

#include <cmath>
#include <cstddef>
#include <set>
#include <vector>

#include "oneapi/tbb.h"
#include "orehov_n_jarvis_pass/common/include/common.hpp"

namespace orehov_n_jarvis_pass {

namespace {

struct BestState {
  Point point;
  bool valid = false;
};

bool IsBetterPoint(Point current, Point candidate, Point best) {
  double orient = OrehovNJarvisPassTBB::CheckLeft(current, best, candidate);
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

BestState FindInitialBestState(const std::vector<Point> &input, Point current) {
  BestState init;
  for (const auto &p : input) {
    if (!(p == current)) {
      init.point = p;
      init.valid = true;
      break;
    }
  }
  return init;
}

BestState ReduceBestStates(const BestState &a, const BestState &b, Point current) {
  if (!a.valid) {
    return b;
  }
  if (!b.valid) {
    return a;
  }
  return BestState{.point = IsBetterPoint(current, b.point, a.point) ? b.point : a.point, .valid = true};
}

}  // namespace

OrehovNJarvisPassTBB::OrehovNJarvisPassTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<Point>();
}

bool OrehovNJarvisPassTBB::ValidationImpl() {
  return (!GetInput().empty());
}

bool OrehovNJarvisPassTBB::PreProcessingImpl() {
  std::set<Point> tmp(GetInput().begin(), GetInput().end());
  input_.assign(tmp.begin(), tmp.end());
  return true;
}

bool OrehovNJarvisPassTBB::RunImpl() {
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

Point OrehovNJarvisPassTBB::FindNext(Point current) const {
  const size_t n = input_.size();

  BestState init = FindInitialBestState(input_, current);
  if (!init.valid) {
    return current;
  }

  auto body = [&](const tbb::blocked_range<size_t> &r, BestState local) -> BestState {
    for (size_t i = r.begin(); i != r.end(); ++i) {
      const Point &p = input_[i];
      if (p == current) {
        continue;
      }
      if (!local.valid || IsBetterPoint(current, p, local.point)) {
        local.point = p;
        local.valid = true;
      }
    }
    return local;
  };

  auto reduce = [&](const BestState &a, const BestState &b) -> BestState { return ReduceBestStates(a, b, current); };

  BestState result = tbb::parallel_reduce(tbb::blocked_range<size_t>(0, n), init, body, reduce);
  return result.point;
}

double OrehovNJarvisPassTBB::CheckLeft(Point a, Point b, Point c) {
  return ((b.x - a.x) * (c.y - a.y)) - ((b.y - a.y) * (c.x - a.x));
}

Point OrehovNJarvisPassTBB::FindFirstElem() const {
  Point current = input_[0];
  for (auto f : input_) {
    if (f.x < current.x || (f.y < current.y && f.x == current.x)) {
      current = f;
    }
  }
  return current;
}

double OrehovNJarvisPassTBB::Distance(Point a, Point b) {
  return std::sqrt(std::pow(a.y - b.y, 2) + std::pow(a.x - b.x, 2));
}

bool OrehovNJarvisPassTBB::PostProcessingImpl() {
  GetOutput() = res_;
  return true;
}

}  // namespace orehov_n_jarvis_pass
