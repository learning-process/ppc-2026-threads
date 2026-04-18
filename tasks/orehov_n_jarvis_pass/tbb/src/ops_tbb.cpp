#include "orehov_n_jarvis_pass/tbb/include/ops_tbb.hpp"

#include <cmath>
#include <cstddef>
#include <set>
#include <vector>

#include "oneapi/tbb.h"
#include "orehov_n_jarvis_pass/common/include/common.hpp"

namespace orehov_n_jarvis_pass {

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
  const auto &input = input_;

  struct Body {
    const Point &current;
    const std::vector<Point> &input;
    Point best_point;

    Body(const Point &c, const std::vector<Point> &in)
        : current(c), input(in), best_point((current == in[0]) ? in[1] : in[0]) {}

    Body(Body &other, tbb::split) : current(other.current), input(other.input), best_point(other.best_point) {}

    void operator()(const tbb::blocked_range<size_t> &range) {
      for (size_t i = range.begin(); i != range.end(); ++i) {
        const Point &point = input[i];
        if (current == point) {
          continue;
        }

        double orient = OrehovNJarvisPassTBB::CheckLeft(current, best_point, point);

        if (orient > 0) {
          best_point = point;
        } else if (orient == 0) {
          if (OrehovNJarvisPassTBB::Distance(current, point) > OrehovNJarvisPassTBB::Distance(current, best_point)) {
            best_point = point;
          }
        }
      }
    }

    void join(const Body &other) {
      double global_orient = OrehovNJarvisPassTBB::CheckLeft(current, best_point, other.best_point);
      if (global_orient > 0) {
        best_point = other.best_point;
      } else if (global_orient == 0) {
        if (OrehovNJarvisPassTBB::Distance(current, other.best_point) >
            OrehovNJarvisPassTBB::Distance(current, best_point)) {
          best_point = other.best_point;
        }
      }
    }
  };

  Body body(current, input);
  tbb::parallel_reduce(tbb::blocked_range<size_t>(0, n), body);

  return body.best_point;
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
