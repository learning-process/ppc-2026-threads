#pragma once

#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace nesterov_a_test_task_threads {

struct Point {
  int x{};
  int y{};
  bool operator==(const Point &other) const = default;
};

struct BinaryImage {
  int width{};
  int height{};
  std::vector<std::uint8_t> data;
};

using InType = BinaryImage;
using OutType = std::vector<std::vector<Point>>;
using TestType = std::tuple<BinaryImage, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace nesterov_a_test_task_threads
