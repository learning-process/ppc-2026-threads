#include "orehov_n_jarvis_pass/omp/include/ops_omp.hpp"

#include <cmath>
#include <set>
#include <vector>

#include "orehov_n_jarvis_pass/common/include/common.hpp"

namespace orehov_n_jarvis_pass {

OrehovNJarvisPassOMP::OrehovNJarvisPassOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<Point>();
}

bool OrehovNJarvisPassOMP::ValidationImpl() {
  return (!GetInput().empty());
}

bool OrehovNJarvisPassOMP::PreProcessingImpl() {
  std::set<Point> tmp(GetInput().begin(), GetInput().end());
  input_.assign(tmp.begin(), tmp.end());
  return true;
}

bool OrehovNJarvisPassOMP::RunImpl() {
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

Point OrehovNJarvisPassOMP::FindLocalBest(Point current, Point initial_next) const {
  Point local_next = initial_next;
  double local_best_orient = -1e9;

  for (const auto &p : input_) {
    if (current == p || local_next == p) {
      continue;
    }

    double orient = CheckLeft(current, local_next, p);

    if (orient > local_best_orient) {
      local_best_orient = orient;
      local_next = p;
    } else if (orient == local_best_orient && orient == 0) {
      if (Distance(current, p) > Distance(current, local_next)) {
        local_next = p;
      }
    }
  }

  return local_next;
}

void OrehovNJarvisPassOMP::UpdateGlobalBest(Point current, Point local_next, Point &global_next) {
  double global_orient = CheckLeft(current, global_next, local_next);
  double current_orient = CheckLeft(current, global_next, global_next);

  if (global_orient > current_orient) {
    global_next = local_next;
  } else if (global_orient == current_orient && global_orient == 0) {
    if (Distance(current, local_next) > Distance(current, global_next)) {
      global_next = local_next;
    }
  }
}

Point OrehovNJarvisPassOMP::FindNext(Point current) const {
  Point next = (current == input_[0]) ? input_[1] : input_[0];

#pragma omp parallel default(none) shared(current, next, input_)
  {
    Point local_next = FindLocalBest(current, next);

#pragma omp critical
    {
      UpdateGlobalBest(current, local_next, next);
    }
  }

  return next;
}

double OrehovNJarvisPassOMP::CheckLeft(Point a, Point b, Point c) {
  return ((b.x - a.x) * (c.y - a.y)) - ((b.y - a.y) * (c.x - a.x));
}

Point OrehovNJarvisPassOMP::FindFirstElem() const {
  Point current = input_[0];
  for (auto f : input_) {
    if (f.x < current.x || (f.y < current.y && f.x == current.x)) {
      current = f;
    }
  }
  return current;
}

double OrehovNJarvisPassOMP::Distance(Point a, Point b) {
  return std::sqrt(pow(a.y - b.y, 2) + pow(a.x - b.x, 2));
}

bool OrehovNJarvisPassOMP::PostProcessingImpl() {
  GetOutput() = res_;
  return true;
}

}  // namespace orehov_n_jarvis_pass
