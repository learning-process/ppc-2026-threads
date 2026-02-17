// tasks/peryashkin_v_binary_component_contour_processing/common/include/common.hpp
#pragma once

#include <cstdint>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace peryashkin_v_binary_component_contour_processing {

struct Point {
  int x{};
  int y{};
};

struct BinaryImage {
  int width{};
  int height{};
  std::vector<std::uint8_t> data;  // 0/1, size = width*height
};

// Вход = бинарное изображение
using InType = BinaryImage;

// Выход = список оболочек по компонентам; каждая оболочка = список вершин
using OutType = std::vector<std::vector<Point>>;

// Для gtest-параметров (если понадобится)
using TestType = std::tuple<int, std::string>;

using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace peryashkin_v_binary_component_contour_processing
