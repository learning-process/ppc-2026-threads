#pragma once

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <ios>
#include <string>
#include <tuple>
#include <vector>

#include "task/include/task.hpp"

namespace ivanova_p_marking_components_on_binary_image {

using InType = int;                // Возвращаем int для совместимости с BaseRunFuncTests
using OutType = std::vector<int>;  // Выходные данные - изображение с метками компонент
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

struct Image {
  int width = 0;
  int height = 0;
  std::vector<uint8_t> data;  // 0 - фон (белый), 1 - объект (черный)
};

// Глобальная переменная для передачи изображения между этапами тестирования
static Image test_image;

}  // namespace ivanova_p_marking_components_on_binary_image
