#pragma once

#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace ilin_a_algorithm_graham {

struct Point {
  double x;
  double y;
};

struct InputData {
  std::vector<Point> points;
};

using InType = InputData;
using OutType = std::vector<Point>;
using TestType = std::tuple<int, InputData, OutType>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace ilin_a_algorithm_graham
