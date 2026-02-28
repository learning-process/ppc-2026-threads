// image_generator.hpp
#pragma once

#include <cstddef>  // Для size_t
#include <cstdint>  // Для uint8_t
#include <fstream>  // Для std::ifstream
#include <ios>      // Для std::ios
#include <string>

#include "ivanova_p_marking_components_on_binary_image/common/include/common.hpp"

namespace ivanova_p_marking_components_on_binary_image {
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

// Вспомогательные функции для каждого тестового случая
namespace test_case_generators {

inline uint8_t GetPixelForTestCase1(int xx, int yy, int width, int height) {
  return test_image_helpers::IsPixelInTestCase1(xx, yy, width, height) ? 1 : 0;
}

inline uint8_t GetPixelForTestCase2(int xx, int yy, int width, int height) {
  return test_image_helpers::IsPixelInTestCase2(xx, yy, width, height) ? 1 : 0;
}

inline uint8_t GetPixelForTestCase3(int xx, int yy, int width, int height) {
  return test_image_helpers::IsPixelInTestCase3(xx, yy, width, height) ? 1 : 0;
}

inline uint8_t GetPixelForTestCase4(int xx, int yy, int width, int height) {
  return test_image_helpers::IsPixelInTestCase4(xx, yy, width, height) ? 1 : 0;
}

inline uint8_t GetPixelForTestCase5(int /*xx*/, int /*yy*/, int /*width*/, int /*height*/) {
  return 0;
}

inline uint8_t GetPixelForTestCase6(int xx, int yy, int width, int height) {
  return (xx == width / 2 && yy == height / 2) ? 1 : 0;
}

inline uint8_t GetPixelForTestCase7(int xx, int yy, int width, int height) {
  return test_image_helpers::IsPixelInTestCase7(xx, yy, width, height) ? 1 : 0;
}

inline uint8_t GetPixelForTestCase8(int xx, int yy, int width, int height) {
  return test_image_helpers::IsPixelInTestCase8(xx, yy, width, height) ? 1 : 0;
}

inline uint8_t GetPixelForTestCase9(int xx, int yy, int width, int height) {
  (void)xx;
  (void)width;
  return (yy == height / 2) ? 1 : 0;
}

inline uint8_t GetPixelForTestCase10(int xx, int yy, int width, int height) {
  (void)yy;
  (void)height;
  return (xx == width / 2) ? 1 : 0;
}

}  // namespace test_case_generators

// Вспомогательная функция для создания тестовых изображений
inline Image CreateTestImage(int width, int height, int test_case) {
  Image img;
  img.width = width;
  img.height = height;
  img.data.resize(static_cast<size_t>(width) * static_cast<size_t>(height));

  // Выбираем функцию-генератор в зависимости от test_case
  uint8_t (*generator)(int, int, int, int) = nullptr;

  switch (test_case) {
    case 1:
      generator = test_case_generators::GetPixelForTestCase1;
      break;
    case 2:
      generator = test_case_generators::GetPixelForTestCase2;
      break;
    case 3:
      generator = test_case_generators::GetPixelForTestCase3;
      break;
    case 4:
      generator = test_case_generators::GetPixelForTestCase4;
      break;
    case 5:
      generator = test_case_generators::GetPixelForTestCase5;
      break;
    case 6:
      generator = test_case_generators::GetPixelForTestCase6;
      break;
    case 7:
      generator = test_case_generators::GetPixelForTestCase7;
      break;
    case 8:
      generator = test_case_generators::GetPixelForTestCase8;
      break;
    case 9:
      generator = test_case_generators::GetPixelForTestCase9;
      break;
    case 10:
      generator = test_case_generators::GetPixelForTestCase10;
      break;
    default:
      generator = test_case_generators::GetPixelForTestCase5;  // фон по умолчанию
      break;
  }

  for (int yy = 0; yy < height; ++yy) {
    for (int xx = 0; xx < width; ++xx) {
      int idx = (yy * width) + xx;
      img.data[idx] = generator(xx, yy, width, height);
    }
  }

  return img;
}
}  // namespace ivanova_p_marking_components_on_binary_image
