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

  struct Body {
    Point current_val;
    const std::vector<Point> *input_ptr;
    Point best_point;

    Body(Point c, const std::vector<Point> *in)
        : current_val(c), input_ptr(in), best_point(GetInitialBestPoint(c, in)) {}

    Body(Body &other, tbb::split /*unused*/)
        : current_val(other.current_val), input_ptr(other.input_ptr), best_point(other.best_point) {}

    static Point GetInitialBestPoint(Point current, const std::vector<Point> *in) {
      if (current == (*in)[0]) {
        return (*in)[1];
      }
      return (*in)[0];
    }

    static bool IsBetter(Point current, Point candidate, Point current_best) {
      double orient = CheckLeft(current, current_best, candidate);
      if (orient > 0) {
        return true;
      }
      if (orient == 0 && Distance(current, candidate) > Distance(current, current_best)) {
        return true;
      }
      return false;
    }

    void operator()(const tbb::blocked_range<size_t> &range) {
      const auto &input = *input_ptr;
      for (size_t i = range.begin(); i != range.end(); ++i) {
        const Point &point = input[i];
        if (current_val == point) {
          continue;
        }
        if (IsBetter(current_val, point, best_point)) {
          best_point = point;
        }
      }
    }

    void Join(const Body &other) {
      if (IsBetter(current_val, other.best_point, best_point)) {
        best_point = other.best_point;
      }
    }

    void join(const Body &other) {
      Join(other);
    }
  };

  Body body(current, &input_);
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
