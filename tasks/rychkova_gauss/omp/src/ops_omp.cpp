#include "rychkova_gauss/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "rychkova_gauss/common/include/common.hpp"
#include "util/include/util.hpp"

namespace rychkova_gauss {

namespace {
int Mirror(int x, int xmin, int xmax) {
  if (x < xmin) {
    return 1;
  }
  if (x >= xmax) {
    return xmax - 1;
  }
  return x;
};

Pixel ComputePixel(const Image &image, std::size_t x, std::size_t y, std::size_t width, std::size_t height) {
  Pixel result = {.R = 0, .G = 0, .B = 0};
  for (int shift_x = -1; shift_x < 2; shift_x++) {
    for (int shift_y = -1; shift_y < 2; shift_y++) {
      int xn = Mirror(static_cast<int>(x) + shift_x, 0, static_cast<int>(width));
      int yn = Mirror(static_cast<int>(y) + shift_y, 0, static_cast<int>(height));
      auto current = image[yn][xn];
      result.R += static_cast<uint8_t>(static_cast<double>(current.R) * kKernel[shift_x + 1][shift_y + 1]);
      result.G += static_cast<uint8_t>(static_cast<double>(current.G) * kKernel[shift_x + 1][shift_y + 1]);
      result.B += static_cast<uint8_t>(static_cast<double>(current.B) * kKernel[shift_x + 1][shift_y + 1]);
    }
  }
  return result;
}
}  // namespace

RychkovaGaussOMP::RychkovaGaussOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = {};  // cписок инициализации - пустой вектор - каждый вложенный вектор и пиксели внутри по умолчанию
}

bool RychkovaGaussOMP::ValidationImpl() {
  if (GetInput().empty()) {
    return false;
  }
  const auto len = GetInput()[0].size();
  return std::ranges::all_of(GetInput(), [len](const auto &row) { return row.size() == len; });
}

bool RychkovaGaussOMP::PreProcessingImpl() {
  return true;
}

bool RychkovaGaussOMP::RunImpl() {
  const auto &image = GetInput();  // сохран.изоб.
  const auto width = image[0].size();
  const auto height = image.size();  // сайз от имаге хранит количествo строк
  GetOutput() = Image(height, std::vector<Pixel>(width, Pixel(0, 0, 0)));
  auto &output = GetOutput();
#pragma omp parallel for shared(width, height, image, output) default(none) num_threads(ppc::util::GetNumThreads())
  for (std::size_t j = 0; j < height; j++) {
    for (std::size_t i = 0; i < width; i++) {
      output[j][i] = ComputePixel(image, i, j, width, height);
    }
  }
  return true;
}

bool RychkovaGaussOMP::PostProcessingImpl() {
  return true;
}

}  // namespace rychkova_gauss
