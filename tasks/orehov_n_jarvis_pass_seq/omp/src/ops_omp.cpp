#include "orehov_n_jarvis_pass_seq/omp/include/ops_omp.hpp"

#include <cmath>
#include <cstddef>
#include <set>
#include <vector>

#include "orehov_n_jarvis_pass_seq/common/include/common.hpp"

namespace orehov_n_jarvis_pass_seq {

OrehovNJarvisPassOMP::OrehovNJarvisPassOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<Point>();
}

bool OrehovNJarvisPassOMP::ValidationImpl() {
  input_ = GetInput();
  return (!GetInput().empty());
}

bool OrehovNJarvisPassOMP::PreProcessingImpl() {
  input_ = GetInput();
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

Point OrehovNJarvisPassOMP::FindNext(Point current) const {
  const auto &input_ref = input_;
  const int n = static_cast<int>(input_ref.size());

  Point initial_candidate = (current == input_ref[0]) ? input_ref[1] : input_ref[0];
  Point global_next = initial_candidate;

#pragma omp parallel default(none) shared(input_ref, n, current, global_next) firstprivate(initial_candidate)
  {
    Point local_next = initial_candidate;

#pragma omp for
    for (int i = 0; i < n; ++i) {
      const Point &point = input_ref[i];
      if (current == point) {
        continue;
      }

      double orient = CheckLeft(current, local_next, point);

      if (orient > 0) {
        local_next = point;
      } else if (std::abs(orient) < 1e-9) {
        if (Distance(current, point) > Distance(current, local_next)) {
          local_next = point;
        }
      }
    }

#pragma omp critical
    {
      double global_orient = CheckLeft(current, global_next, local_next);
      if (global_orient > 0) {
        global_next = local_next;
      } else if (std::abs(global_orient) < 1e-9) {
        if (Distance(current, local_next) > Distance(current, global_next)) {
          global_next = local_next;
        }
      }
    }
  }

  return global_next;
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

}  // namespace orehov_n_jarvis_pass_seq
