#include "romanov_a_gauss_block/seq/include/ops_seq.hpp"

#include <numeric>
#include <vector>
#include <cstdint>
#include <tuple>

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

bool RomanovAGaussBlockSEQ::RunImpl() {

  const int width = std::get<0>(GetInput());
  const int height = std::get<1>(GetInput());

  const std::vector<uint8_t> initial_picture = std::get<2>(GetInput());
  std::vector<uint8_t> result_picture(height * width * 3);

  const int kernel[3][3] = {
    {1, 2, 1},
    {2, 4, 2},
    {1, 2, 1}
  };

  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      for (int c = 0; c < 3; ++c) {
        int sum = 0;
        for (int dy = -1; dy <= 1; ++dy) {
          for (int dx = -1; dx <= 1; ++dx) {
            int nx = x + dx;
            int ny = y + dy;
            if (nx < 0) {
              nx = 0;
            }
            if (nx >= width) {
              nx = width - 1;
            }
            if (ny < 0) {
              ny = 0;
            }
            if (ny >= height) {
              ny = height - 1;
            }

            uint8_t pixel = initial_picture[(((ny * width) + nx)) * 3 + c];
            sum += static_cast<int>(pixel) * kernel[dy + 1][dx + 1];
          }
        }

        int result_value = (sum + 8) / 16;
        if (result_value < 0) {
          result_value = 0;
        }
        if (result_value > 255) {
          result_value = 255;
        }
        result_picture[(((y * width) + x) * 3) + c] = static_cast<uint8_t>(result_value);
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
