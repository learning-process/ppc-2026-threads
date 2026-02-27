#include "kutuzov_i_convex_hull_jarvis/seq/include/ops_seq.hpp"

#include <cmath>
#include <vector>

#include "kutuzov_i_convex_hull_jarvis/common/include/common.hpp"

namespace kutuzov_i_convex_hull_jarvis {

KutuzovITestConvexHullSEQ::KutuzovITestConvexHullSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};
}

double KutuzovITestConvexHullSEQ::DistanceSquared(double a_x, double a_y, double b_x, double b_y) {
  return ((a_x - b_x) * (a_x - b_x)) + ((a_y - b_y) * (a_y - b_y));
}

double KutuzovITestConvexHullSEQ::CrossProduct(double o_x, double o_y, double a_x, double a_y, double b_x, double b_y) {
  return ((a_x - o_x) * (b_y - o_y)) - ((a_y - o_y) * (b_x - o_x));
}

bool KutuzovITestConvexHullSEQ::IsBetterPoint(double cross, double epsilon, double current_x, double current_y,
                                              double i_x, double i_y, double next_x, double next_y) {
  if (cross < -epsilon) {
    return true;
  }

  if (std::abs(cross) < epsilon) {
    return DistanceSquared(current_x, current_y, i_x, i_y) > DistanceSquared(current_x, current_y, next_x, next_y);
  }

  return false;
}

bool KutuzovITestConvexHullSEQ::ValidationImpl() {
  return true;
}

bool KutuzovITestConvexHullSEQ::PreProcessingImpl() {
  return true;
}

bool KutuzovITestConvexHullSEQ::RunImpl() {
  if (GetInput().size() < 3) {
    GetOutput() = GetInput();
    return true;
  }

  // Finding left-most point
  int leftmost = 0;
  double leftmost_x = std::get<0>(GetInput()[leftmost]);
  double leftmost_y = std::get<1>(GetInput()[leftmost]);

  for (size_t i = 0; i < GetInput().size(); ++i) {
    double x = std::get<0>(GetInput()[i]);
    double y = std::get<1>(GetInput()[i]);

    if ((x < leftmost_x) || ((x == leftmost_x) && (y < leftmost_y))) {
      leftmost = i;
      leftmost_x = std::get<0>(GetInput()[leftmost]);
      leftmost_y = std::get<1>(GetInput()[leftmost]);
    }
  }

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

    for (size_t i = 0; i < GetInput().size(); ++i) {
      if (i == current) {
        continue;
      }

      double i_x = std::get<0>(GetInput()[i]);
      double i_y = std::get<1>(GetInput()[i]);

      double cross = CrossProduct(current_x, current_y, next_x, next_y, i_x, i_y);

      if (IsBetterPoint(cross, epsilon, current_x, current_y, i_x, i_y, next_x, next_y)) {
        next = i;
        next_x = std::get<0>(GetInput()[next]);
        next_y = std::get<1>(GetInput()[next]);
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

bool KutuzovITestConvexHullSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace kutuzov_i_convex_hull_jarvis
