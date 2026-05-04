#include "kutuzov_i_convex_hull_jarvis/omp/include/ops_omp.hpp"

#include <omp.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include "kutuzov_i_convex_hull_jarvis/common/include/common.hpp"

namespace kutuzov_i_convex_hull_jarvis {

KutuzovITestConvexHullOMP::KutuzovITestConvexHullOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

double KutuzovITestConvexHullOMP::DistanceSquared(double a_x, double a_y, double b_x, double b_y) {
  return ((a_x - b_x) * (a_x - b_x)) + ((a_y - b_y) * (a_y - b_y));
}

double KutuzovITestConvexHullOMP::CrossProduct(double o_x, double o_y, double a_x, double a_y, double b_x, double b_y) {
  return ((a_x - o_x) * (b_y - o_y)) - ((a_y - o_y) * (b_x - o_x));
}

size_t KutuzovITestConvexHullOMP::FindLeftmostPoint(const InType &input) {
  size_t leftmost = 0;

#pragma omp parallel default(none) shared(input, leftmost)
  {
    size_t local_leftmost = 0;
    double local_leftmost_x = std::get<0>(input[local_leftmost]);
    double local_leftmost_y = std::get<1>(input[local_leftmost]);

#pragma omp for nowait
    for (size_t i = 0; i < input.size(); ++i) {
      double x = std::get<0>(input[i]);
      double y = std::get<1>(input[i]);

      if ((x < local_leftmost_x) || ((x == local_leftmost_x) && (y < local_leftmost_y))) {
        local_leftmost = i;
        local_leftmost_x = x;
        local_leftmost_y = y;
      }
    }

#pragma omp critical
    {
      double leftmost_x = std::get<0>(input[leftmost]);
      double leftmost_y = std::get<1>(input[leftmost]);

      if ((local_leftmost_x < leftmost_x) || ((local_leftmost_x == leftmost_x) && (local_leftmost_y < leftmost_y))) {
        leftmost = local_leftmost;
      }
    }
  }
  return leftmost;
}

bool KutuzovITestConvexHullOMP::IsBetterPoint(double cross, double epsilon, double current_x, double current_y,
                                              double i_x, double i_y, double next_x, double next_y) {
  if (cross < -epsilon) {
    return true;
  }

  if (std::abs(cross) < epsilon) {
    return DistanceSquared(current_x, current_y, i_x, i_y) > DistanceSquared(current_x, current_y, next_x, next_y);
  }

  return false;
}

bool KutuzovITestConvexHullOMP::ValidationImpl() {
  return true;
}

bool KutuzovITestConvexHullOMP::PreProcessingImpl() {
  return true;
}

bool KutuzovITestConvexHullOMP::RunImpl() {
  if (GetInput().size() < 3) {
    GetOutput() = GetInput();
    return true;
  }

  // Finding left-most point
  size_t leftmost = FindLeftmostPoint(GetInput());

  // Main loop
  size_t current = leftmost;
  double current_x = std::get<0>(GetInput()[current]);
  double current_y = std::get<1>(GetInput()[current]);

  const double epsilon = 1e-9;

  while (true) {
    // Adding current point to the hull
    GetOutput().push_back(GetInput()[current]);

    // Finding the next point of the hull
    size_t next = (current + 1) % GetInput().size();
    double next_x = std::get<0>(GetInput()[next]);
    double next_y = std::get<1>(GetInput()[next]);

#pragma omp parallel default(none) shared(current, current_x, current_y, next, next_x, next_y, epsilon)
    {
      size_t local_next = (current + 1) % GetInput().size();
      double local_next_x = std::get<0>(GetInput()[local_next]);
      double local_next_y = std::get<1>(GetInput()[local_next]);

#pragma omp for nowait
      for (size_t i = 0; i < GetInput().size(); ++i) {
        if (i == current) {
          continue;
        }

        double i_x = std::get<0>(GetInput()[i]);
        double i_y = std::get<1>(GetInput()[i]);

        double cross = CrossProduct(current_x, current_y, local_next_x, local_next_y, i_x, i_y);

        if (IsBetterPoint(cross, epsilon, current_x, current_y, i_x, i_y, local_next_x, local_next_y)) {
          local_next = i;
          local_next_x = i_x;
          local_next_y = i_y;
        }
      }

#pragma omp critical
      {
        double cross = CrossProduct(current_x, current_y, next_x, next_y, local_next_x, local_next_y);

        if (IsBetterPoint(cross, epsilon, current_x, current_y, local_next_x, local_next_y, next_x, next_y)) {
          next = local_next;
          next_x = local_next_x;
          next_y = local_next_y;
        }
      }
    }

    current = next;
    current_x = next_x;
    current_y = next_y;

    // Loop until we wrap around to the first point
    if (current == leftmost) {
      break;
    }
  }
  return true;
}

bool KutuzovITestConvexHullOMP::PostProcessingImpl() {
  return true;
}

}  // namespace kutuzov_i_convex_hull_jarvis
