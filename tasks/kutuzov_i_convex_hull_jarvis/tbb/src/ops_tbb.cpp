#include "kutuzov_i_convex_hull_jarvis/tbb/include/ops_tbb.hpp"

#include <tbb/blocked_range.h>
#include <tbb/combinable.h>
#include <tbb/parallel_for.h>
#include <tbb/partitioner.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include "kutuzov_i_convex_hull_jarvis/common/include/common.hpp"

namespace kutuzov_i_convex_hull_jarvis {

namespace {

struct BestPoint {
  size_t idx = 0;
  double x = 0.0;
  double y = 0.0;
};

size_t FindLeftmostPointHelper(const InType &input) {
  BestPoint init;
  init.idx = 0;
  init.x = std::get<0>(input[0]);
  init.y = std::get<1>(input[0]);

  tbb::combinable<BestPoint> local_best([&]() { return init; });

  tbb::parallel_for(tbb::blocked_range<size_t>(0, input.size()), [&](const tbb::blocked_range<size_t> &r) {
    BestPoint &best = local_best.local();
    for (size_t i = r.begin(); i != r.end(); ++i) {
      double ix = std::get<0>(input[i]);
      double iy = std::get<1>(input[i]);
      if ((ix < best.x) || ((ix == best.x) && (iy < best.y))) {
        best.idx = i;
        best.x = ix;
        best.y = iy;
      }
    }
  }, tbb::static_partitioner{});

  BestPoint global = local_best.combine([](const BestPoint &a, const BestPoint &b) -> BestPoint {
    if ((b.x < a.x) || ((b.x == a.x) && (b.y < a.y))) {
      return b;
    }
    return a;
  });
  return global.idx;
}

size_t FindNextPointHelper(const InType &input, size_t current, double current_x, double current_y, double epsilon) {
  size_t init_idx = (current + 1) % input.size();
  BestPoint init;
  init.idx = init_idx;
  init.x = std::get<0>(input[init_idx]);
  init.y = std::get<1>(input[init_idx]);

  tbb::combinable<BestPoint> local_next([&]() { return init; });

  tbb::parallel_for(tbb::blocked_range<size_t>(0, input.size()), [&](const tbb::blocked_range<size_t> &r) {
    BestPoint &best = local_next.local();
    for (size_t i = r.begin(); i != r.end(); ++i) {
      if (i == current) {
        continue;
      }

      double ix = std::get<0>(input[i]);
      double iy = std::get<1>(input[i]);

      double cross = KutuzovITestConvexHullTBB::CrossProduct(current_x, current_y, best.x, best.y, ix, iy);

      if (KutuzovITestConvexHullTBB::IsBetterPoint(cross, epsilon, current_x, current_y, ix, iy, best.x, best.y)) {
        best.idx = i;
        best.x = ix;
        best.y = iy;
      }
    }
  }, tbb::static_partitioner{});

  BestPoint global = local_next.combine([&](const BestPoint &a, const BestPoint &b) -> BestPoint {
    double cross = KutuzovITestConvexHullTBB::CrossProduct(current_x, current_y, a.x, a.y, b.x, b.y);
    if (KutuzovITestConvexHullTBB::IsBetterPoint(cross, epsilon, current_x, current_y, b.x, b.y, a.x, a.y)) {
      return b;
    }

    return a;
  });
  return global.idx;
}

}  // anonymous namespace

KutuzovITestConvexHullTBB::KutuzovITestConvexHullTBB(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

double KutuzovITestConvexHullTBB::DistanceSquared(double a_x, double a_y, double b_x, double b_y) {
  return ((a_x - b_x) * (a_x - b_x)) + ((a_y - b_y) * (a_y - b_y));
}

double KutuzovITestConvexHullTBB::CrossProduct(double o_x, double o_y, double a_x, double a_y, double b_x, double b_y) {
  return ((a_x - o_x) * (b_y - o_y)) - ((a_y - o_y) * (b_x - o_x));
}

bool KutuzovITestConvexHullTBB::IsBetterPoint(double cross, double epsilon, double current_x, double current_y,
                                              double i_x, double i_y, double next_x, double next_y) {
  if (cross < -epsilon) {
    return true;
  }
  if (std::abs(cross) < epsilon) {
    return DistanceSquared(current_x, current_y, i_x, i_y) > DistanceSquared(current_x, current_y, next_x, next_y);
  }
  return false;
}

bool KutuzovITestConvexHullTBB::ValidationImpl() {
  return true;
}

bool KutuzovITestConvexHullTBB::PreProcessingImpl() {
  return true;
}

bool KutuzovITestConvexHullTBB::RunImpl() {
  if (GetInput().size() < 3) {
    GetOutput() = GetInput();
    return true;
  }

  size_t leftmost = FindLeftmostPointHelper(GetInput());
  size_t current = leftmost;
  double current_x = std::get<0>(GetInput()[current]);
  double current_y = std::get<1>(GetInput()[current]);
  const double epsilon = 1e-9;

  while (current != leftmost || GetOutput().empty()) {
    GetOutput().push_back(GetInput()[current]);

    size_t next = FindNextPointHelper(GetInput(), current, current_x, current_y, epsilon);

    current = next;
    current_x = std::get<0>(GetInput()[current]);
    current_y = std::get<1>(GetInput()[current]);
  }
  return true;
}

bool KutuzovITestConvexHullTBB::PostProcessingImpl() {
  return true;
}

}  // namespace kutuzov_i_convex_hull_jarvis
