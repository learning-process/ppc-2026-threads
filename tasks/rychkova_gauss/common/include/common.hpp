#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "task/include/task.hpp"

namespace rychkova_gauss {

struct Pixel {
  uint8_t R;
  uint8_t G;
  uint8_t B;
  bool operator==(const Pixel &other) const {
    return R == other.R && G == other.G && B == other.B;
  }
};
const int kEPS = 1;

const std::vector<std::vector<double>> kKernel = {
    {1. / 16, 2. / 16, 1. / 16}, {2. / 16, 4. / 16, 2. / 16}, {1. / 16, 2. / 16, 1. / 16}};
using Image = std::vector<std::vector<Pixel>>;  // создали псевдоним для изображения

using InType = Image;          // вход
using OutType = Image;         // выход
using TestType = std::string;  // tuрle нужен только потому что работает с гугл
using BaseTask = ppc::task::Task<InType, OutType>;

}  // namespace rychkova_gauss
