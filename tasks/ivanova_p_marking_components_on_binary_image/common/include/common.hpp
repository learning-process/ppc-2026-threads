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

// Функция для загрузки изображения из текстового файла
inline Image LoadImageFromTxt(const std::string &filename) {
  Image img;
  std::ifstream file(filename);

  if (!file.is_open()) {
    img.width = 0;
    img.height = 0;
    return img;
  }

  file >> img.width >> img.height;

  if (img.width <= 0 || img.height <= 0) {
    img.width = 0;
    img.height = 0;
    return img;
  }

  img.data.resize(static_cast<size_t>(img.width) * static_cast<size_t>(img.height));

  for (size_t i = 0; i < img.data.size(); ++i) {
    int pixel_value = 0;
    if (!(file >> pixel_value)) {
      img.width = 0;
      img.height = 0;
      return img;
    }
    // Преобразуем в бинарное: если пиксель темный (< 128), то 1, иначе 0
    // 0 - объект (черный), 255 - фон (белый)
    img.data[i] = (pixel_value < 128) ? 1 : 0;
  }

  file.close();
  return img;
}

// Функция для загрузки PPM изображения
inline Image LoadPPMImage(const std::string &filename) {
  Image img;
  std::ifstream file(filename, std::ios::binary);

  if (!file.is_open()) {
    img.width = 0;
    img.height = 0;
    return img;
  }

  std::string magic;
  file >> magic;

  if (magic != "P5" && magic != "P6") {
    img.width = 0;
    img.height = 0;
    return img;
  }

  int max_val = 0;
  file >> img.width >> img.height >> max_val;
  file.get();  // Пропускаем пробел после max_val

  img.data.resize(static_cast<size_t>(img.width) * static_cast<size_t>(img.height));

  if (magic == "P5") {
    // Grayscale
    std::vector<uint8_t> gray_data(static_cast<size_t>(img.width) * static_cast<size_t>(img.height));
    file.read(reinterpret_cast<char *>(gray_data.data()),
              static_cast<std::streamsize>(static_cast<size_t>(img.width) * static_cast<size_t>(img.height)));

    for (size_t i = 0; i < gray_data.size(); ++i) {
      // Преобразуем в бинарное: если пиксель темный (< 128), то 1, иначе 0
      img.data[i] = (gray_data[i] < 128) ? 1 : 0;
    }
  } else {
    // RGB
    std::vector<uint8_t> rgb_data(static_cast<size_t>(img.width) * static_cast<size_t>(img.height) * 3U);
    file.read(reinterpret_cast<char *>(rgb_data.data()),
              static_cast<std::streamsize>(static_cast<size_t>(img.width) * static_cast<size_t>(img.height) * 3U));

    for (size_t i = 0; i < img.data.size(); ++i) {
      uint8_t r = rgb_data[i * 3U];
      uint8_t g = rgb_data[(i * 3U) + 1U];
      uint8_t b = rgb_data[(i * 3U) + 2U];

      // Преобразуем в бинарное: если пиксель темный, то 1, иначе 0
      uint8_t gray = (r + g + b) / 3;
      img.data[i] = (gray < 128) ? 1 : 0;
    }
  }

  file.close();
  return img;
}

// Вспомогательные функции для создания тестовых изображений
namespace test_image_helpers {

inline bool IsPixelInTestCase1(int xx, int yy, int width, int height) {
  return xx > width / 4 && xx < (3 * width) / 4 && yy > height / 4 && yy < (3 * height) / 4;
}

inline bool IsPixelInTestCase2(int xx, int yy, int width, int height) {
  return (xx > width / 8 && xx < (3 * width) / 8 && yy > height / 8 && yy < (3 * height) / 8) ||
         (xx > (5 * width) / 8 && xx < (7 * width) / 8 && yy > (5 * height) / 8 && yy < (7 * height) / 8);
}

inline bool IsPixelInTestCase3(int xx, int yy, int width, int height) {
  return (xx > width / 10 && xx < (3 * width) / 10 && yy > height / 10 && yy < (3 * height) / 10) ||
         (xx > (4 * width) / 10 && xx < (6 * width) / 10 && yy > (4 * height) / 10 && yy < (6 * height) / 10) ||
         (xx > (7 * width) / 10 && xx < (9 * width) / 10 && yy > (7 * height) / 10 && yy < (9 * height) / 10);
}

inline bool IsPixelInTestCase4(int xx, int yy, int width, int height) {
  return (xx > width / 3 && xx < ((width / 3) + 5) && yy > height / 4 && yy < (3 * height) / 4) ||
         (xx > (2 * width) / 3 && xx < (((2 * width) / 3) + 5) && yy > height / 4 && yy < (3 * height) / 4) ||
         (xx > width / 3 && xx < (((2 * width) / 3) + 5) && yy > ((height / 2) - 2) && yy < ((height / 2) + 2));
}

inline bool IsPixelInTestCase7(int xx, int yy, int width, int height) {
  return (xx == width / 2 && yy == height / 4) || (xx == (3 * width) / 4 && yy == height / 4) ||
         (xx == width / 4 && yy == (3 * height) / 4) || (xx == (3 * width) / 4 && yy == (3 * height) / 4);
}

inline bool IsPixelInTestCase8(int xx, int yy, int width, int height) {
  int cell_width = width / 3;
  int cell_height = height / 3;
  int local_x = xx % cell_width;
  int local_y = yy % cell_height;
  return local_x > cell_width / 4 && local_x < (3 * cell_width) / 4 && local_y > cell_height / 4 &&
         local_y < (3 * cell_height) / 4;
}

}  // namespace test_image_helpers

// Вспомогательная функция для создания тестовых изображений
inline Image CreateTestImage(int width, int height, int test_case) {
  Image img;
  img.width = width;
  img.height = height;
  img.data.resize(static_cast<size_t>(width) * static_cast<size_t>(height));

  for (int yy = 0; yy < height; ++yy) {
    for (int xx = 0; xx < width; ++xx) {
      int idx = (yy * width) + xx;
      uint8_t pixel = 0;  // фон по умолчанию

      switch (test_case) {
        case 1:
          pixel = test_image_helpers::IsPixelInTestCase1(xx, yy, width, height) ? 1 : 0;
          break;
        case 2:
          pixel = test_image_helpers::IsPixelInTestCase2(xx, yy, width, height) ? 1 : 0;
          break;
        case 3:
          pixel = test_image_helpers::IsPixelInTestCase3(xx, yy, width, height) ? 1 : 0;
          break;
        case 4:
          pixel = test_image_helpers::IsPixelInTestCase4(xx, yy, width, height) ? 1 : 0;
          break;
        case 5:
          pixel = 0;
          break;
        case 6:
          pixel = (xx == width / 2 && yy == height / 2) ? 1 : 0;
          break;
        case 7:
          pixel = test_image_helpers::IsPixelInTestCase7(xx, yy, width, height) ? 1 : 0;
          break;
        case 8:
          pixel = test_image_helpers::IsPixelInTestCase8(xx, yy, width, height) ? 1 : 0;
          break;
        case 9:
          pixel = (yy == height / 2) ? 1 : 0;
          break;
        case 10:
          pixel = (xx == width / 2) ? 1 : 0;
          break;
        default:
          pixel = 0;
          break;
      }

      img.data[idx] = pixel;
    }
  }

  return img;
}

}  // namespace ivanova_p_marking_components_on_binary_image
