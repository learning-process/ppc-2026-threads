#include "belov_e_sobel/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "belov_e_sobel/common/include/common.hpp"

namespace belov_e_sobel {

BelovESobelSEQ::BelovESobelSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = in;
}

bool BelovESobelSEQ::ValidationImpl() {
  return !std::get<0>(GetInput()).empty() && (std::get<1>(GetInput()) > 0) && (std::get<2>(GetInput()) > 0);
}

bool BelovESobelSEQ::PreProcessingImpl() {
  return true;
}

bool BelovESobelSEQ::RunImpl() {
  const std::vector<uint8_t> &input = std::get<0>(GetInput());
  std::vector<uint8_t> &output = std::get<0>(GetOutput());
  int width = std::get<1>(GetInput());
  int height = std::get<2>(GetInput());

  auto get_px = [&](int col, int row) -> float {
    int clamped_x = std::clamp(col, 0, width - 1);
    int clamped_y = std::clamp(row, 0, height - 1);
    return static_cast<float>(input[(static_cast<size_t>(clamped_y) * width) + clamped_x]);
  };

  for (int row = 0; row < height; ++row) {
    for (int col = 0; col < width; ++col) {
      // ядро Gx (горизонтальное)
      float gx = (-1 * get_px(col - 1, row - 1)) + (1 * get_px(col + 1, row - 1)) + (-2 * get_px(col - 1, row)) +
                 (2 * get_px(col + 1, row)) + (-1 * get_px(col - 1, row + 1)) + (1 * get_px(col + 1, row + 1));
      // ядро Gy (вертикальное)
      float gy = (-1 * get_px(col - 1, row - 1)) - (2 * get_px(col, row - 1)) - (1 * get_px(col + 1, row - 1)) +
                 (1 * get_px(col - 1, row + 1)) + (2 * get_px(col, row + 1)) + (1 * get_px(col + 1, row + 1));
      // Результат
      float magnitude = std::sqrt((gx * gx) + (gy * gy));
      output[(static_cast<size_t>(row) * width) + col] = static_cast<uint8_t>(std::min(255.0F, magnitude));
    }
  }

  return true;
}

bool BelovESobelSEQ::PostProcessingImpl() {
  return !std::get<0>(GetOutput()).empty() && (std::get<1>(GetOutput()) > 0) && (std::get<2>(GetOutput()) > 0);
}

}  // namespace belov_e_sobel
