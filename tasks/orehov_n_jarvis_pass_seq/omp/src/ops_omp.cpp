#include "orehov_n_jarvis_pass_seq/omp/include/ops_omp.hpp"

#include <cmath>
#include <set>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "orehov_n_jarvis_pass_seq/common/include/common.hpp"

namespace orehov_n_jarvis_pass_seq {

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
    Point next = FindNextOMP(current);
    if (next == res_[0]) {
      break;
    }

    current = next;
    res_.push_back(next);
  }

  return true;
}

Point OrehovNJarvisPassOMP::FindNextOMP(Point current) const {
  Point next = current == input_[0] ? input_[1] : input_[0];
  std::vector<Point> input_copy = input_;
  
  #pragma omp parallel default(none) shared(current, input_copy, next)
  {
    Point local_next = next;
    
    #pragma omp for nowait
    for (const auto &p : input_copy) {
      if (current == p || local_next == p) {
        continue;
      }
      
      double orient = CheckLeft(current, local_next, p);
      if (orient > 0) {
        local_next = p;
      } else if (orient == 0) {
        if (Distance(current, local_next) < Distance(current, p)) {
          local_next = p;
        }
      }
    }
    
    #pragma omp critical
    {
      double orient = CheckLeft(current, next, local_next);
      if (orient > 0) {
        next = local_next;
      } else if (orient == 0) {
        if (Distance(current, next) < Distance(current, local_next)) {
          next = local_next;
        }
      }
    }
  }
  
  return next;
}

double OrehovNJarvisPassOMP::CheckLeft(Point a, Point b, Point c) {
  return ((b.x - a.x) * (c.y - a.y)) - ((b.y - a.y) * (c.x - a.x));
}

Point OrehovNJarvisPassOMP::FindFirstElem() const {
  Point current = input_[0];
  std::vector<Point> input_copy = input_;
  
  #pragma omp parallel default(none) shared(input_copy, current)
  {
    Point local_min = current;
    
    #pragma omp for nowait
    for (const auto &f : input_copy) {
      if (f.x < local_min.x || (f.y < local_min.y && f.x == local_min.x)) {
        local_min = f;
      }
    }
    
    #pragma omp critical
    {
      if (local_min.x < current.x || (local_min.y < current.y && local_min.x == current.x)) {
        current = local_min;
      }
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

}  // namespace orehov_n_jarvis_pass_seq