#pragma once

#include <string>
#include <tuple>
#include <vector>
#include <cstdint>
#include <fstream>
#include <sstream>

#include "task/include/task.hpp"

namespace ivanova_p_marking_components_on_binary_image {

using InType = int;  // Возвращаем int для совместимости с BaseRunFuncTests
using OutType = std::vector<int>;  // Выходные данные - изображение с метками компонент
using TestType = std::tuple<int, std::string>;
using BaseTask = ppc::task::Task<InType, OutType>;

struct Image {
    int width;
    int height;
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
  
  img.data.resize(img.width * img.height);
  
  for (int i = 0; i < img.width * img.height; ++i) {
    int pixel_value;
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
  
  int max_val;
  file >> img.width >> img.height >> max_val;
  file.get();  // Пропускаем пробел после max_val
  
  img.data.resize(img.width * img.height);
  
  if (magic == "P5") {
    // Grayscale
    std::vector<uint8_t> gray_data(img.width * img.height);
    file.read(reinterpret_cast<char*>(gray_data.data()), img.width * img.height);
    
    for (int i = 0; i < img.width * img.height; ++i) {
      // Преобразуем в бинарное: если пиксель темный (< 128), то 1, иначе 0
      img.data[i] = (gray_data[i] < 128) ? 1 : 0;
    }
  } else {
    // RGB
    std::vector<uint8_t> rgb_data(img.width * img.height * 3);
    file.read(reinterpret_cast<char*>(rgb_data.data()), img.width * img.height * 3);
    
    for (int i = 0; i < img.width * img.height; ++i) {
      uint8_t r = rgb_data[i * 3];
      uint8_t g = rgb_data[i * 3 + 1];
      uint8_t b = rgb_data[i * 3 + 2];
      
      // Преобразуем в бинарное: если пиксель темный, то 1, иначе 0
      uint8_t gray = (r + g + b) / 3;
      img.data[i] = (gray < 128) ? 1 : 0;
    }
  }
  
  file.close();
  return img;
}

// Вспомогательная функция для создания тестовых изображений
inline Image CreateTestImage(int width, int height, int test_case) {
  Image img;
  img.width = width;
  img.height = height;
  img.data.resize(width * height);
  
  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      int idx = y * width + x;
      uint8_t pixel = 0;  // фон по умолчанию
      
      switch (test_case) {
        case 1: // Один прямоугольник в центре
          if (x > width/4 && x < 3*width/4 && y > height/4 && y < 3*height/4) {
            pixel = 1;
          }
          break;
          
        case 2: // Два отдельных прямоугольника
          if ((x > width/8 && x < 3*width/8 && y > height/8 && y < 3*height/8) ||
              (x > 5*width/8 && x < 7*width/8 && y > 5*height/8 && y < 7*height/8)) {
            pixel = 1;
          }
          break;
          
        case 3: // Три отдельных компонента
          if ((x > width/10 && x < 3*width/10 && y > height/10 && y < 3*height/10) ||
              (x > 4*width/10 && x < 6*width/10 && y > 4*height/10 && y < 6*height/10) ||
              (x > 7*width/10 && x < 9*width/10 && y > 7*height/10 && y < 9*height/10)) {
            pixel = 1;
          }
          break;
          
        case 4: // Связанные компоненты (буква "H")
          if ((x > width/3 && x < width/3 + 5 && y > height/4 && y < 3*height/4) ||
              (x > 2*width/3 && x < 2*width/3 + 5 && y > height/4 && y < 3*height/4) ||
              (x > width/3 && x < 2*width/3 + 5 && y > height/2 - 2 && y < height/2 + 2)) {
            pixel = 1;
          }
          break;
          
        case 5: // Все фон (нет объектов)
          pixel = 0;
          break;
          
        case 6: // Одиночный пиксель в центре
          if (x == width/2 && y == height/2) {
            pixel = 1;
          }
          break;
          
        case 7: // Диагональные соседи (не должны объединяться)
          // 4 пикселя по диагонали - каждый отдельная компонента
          if ((x == width/4 && y == height/4) ||
              (x == 3*width/4 && y == height/4) ||
              (x == width/4 && y == 3*height/4) ||
              (x == 3*width/4 && y == 3*height/4)) {
            pixel = 1;
          }
          break;
          
        case 8: // 9 маленьких компонент (3x3 сетка)
          {
            int cell_width = width / 3;
            int cell_height = height / 3;
            int local_x = x % cell_width;
            int local_y = y % cell_height;
            
            // Маленький квадрат в центре каждой ячейки
            if (local_x > cell_width/4 && local_x < 3*cell_width/4 &&
                local_y > cell_height/4 && local_y < 3*cell_height/4) {
              pixel = 1;
            }
          }
          break;
          
        case 9: // Горизонтальная линия в центре
          if (y == height/2) {
            pixel = 1;
          }
          break;
          
        case 10: // Вертикальная линия в центре
          if (x == width/2) {
            pixel = 1;
          }
          break;
      }
      
      img.data[idx] = pixel;
    }
  }
  
  return img;
}

}  // namespace ivanova_p_marking_components_on_binary_image