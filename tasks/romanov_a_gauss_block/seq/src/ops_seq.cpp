#include "romanov_a_gauss_block/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstdint>
#include <tuple>
#include <vector>

#include "romanov_a_gauss_block/common/include/common.hpp"

namespace romanov_a_gauss_block {

RomanovAGaussBlockSEQ::RomanovAGaussBlockSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<uint8_t>();
}

bool RomanovAGaussBlockSEQ::ValidationImpl() {
  return std::get<0>(GetInput()) * std::get<1>(GetInput()) * 3 == static_cast<int>(std::get<2>(GetInput()).size());
}

bool RomanovAGaussBlockSEQ::PreProcessingImpl() {
  return true;
}

namespace {
int ApplyKernel(const std::vector<uint8_t> &img, int row, int col, int channel, int width, int height,
                const std::array<std::array<int, 3>, 3> &kernel) {
  int sum = 0;
  for (int dr = -1; dr <= 1; ++dr) {
    for (int dc = -1; dc <= 1; ++dc) {
      int nr = row + dr;
      int nc = col + dc;
      if (nr >= 0 && nr < height && nc >= 0 && nc < width) {
        size_t idx = (static_cast<size_t>(((nr * width) + nc) * 3) + channel);
        sum += static_cast<int>(img[idx]) * kernel[dr + 1][dc + 1];
      }
    }
  }
  return sum;
}
}  // namespace

bool RomanovAGaussBlockSEQ::RunImpl() {
  const int width = std::get<0>(GetInput());
  const int height = std::get<1>(GetInput());

  const std::vector<uint8_t> initial_picture = std::get<2>(GetInput());
  std::vector<uint8_t> result_picture(height * width * 3);

  const std::array<std::array<int, 3>, 3> kernel = {{{1, 2, 1}, {2, 4, 2}, {1, 2, 1}}};

  for (int row = 0; row < height; ++row) {
    for (int col = 0; col < width; ++col) {
      for (int channel = 0; channel < 3; ++channel) {
        int sum = ApplyKernel(initial_picture, row, col, channel, width, height, kernel);
        int result_value = (sum + 8) / 16;
        result_value = std::clamp(result_value, 0, 255);
        size_t idx = static_cast<size_t>((((row * width) + col) * 3) + channel);
        result_picture[idx] = static_cast<uint8_t>(result_value);
      }
    }
  }

  GetOutput() = result_picture;
  return true;
}

bool RomanovAGaussBlockSEQ::PostProcessingImpl() {
  return true;
}

}  // namespace romanov_a_gauss_block
