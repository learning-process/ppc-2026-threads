#include "romanov_a_gauss_block/seq/include/ops_seq.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
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
  for (size_t kr = 0; kr < 3; ++kr) {
    for (size_t kc = 0; kc < 3; ++kc) {
      int nr = row + static_cast<int>(kr) - 1;
      int nc = col + static_cast<int>(kc) - 1;
      if (nr >= 0 && nr < height && nc >= 0 && nc < width) {
        size_t idx = (((static_cast<size_t>(nr) * width) + nc) * 3) + channel;
        sum += (static_cast<int>(img[idx]) * kernel.at(kr).at(kc));
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
  std::vector<uint8_t> result_picture(static_cast<size_t>(height * width * 3));

  const std::array<std::array<int, 3>, 3> kernel = {{{1, 2, 1}, {2, 4, 2}, {1, 2, 1}}};

  for (int row = 0; row < height; ++row) {
    for (int col = 0; col < width; ++col) {
      for (int channel = 0; channel < 3; ++channel) {
        int sum = ApplyKernel(initial_picture, row, col, channel, width, height, kernel);
        int result_value = (sum + 8) / 16;
        result_value = std::clamp(result_value, 0, 255);
        auto idx = ((static_cast<size_t>(row) * width + col) * 3) + channel;
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
