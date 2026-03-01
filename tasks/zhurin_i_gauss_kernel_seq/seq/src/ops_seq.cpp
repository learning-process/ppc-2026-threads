#include "zhurin_i_gauss_kernel_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

namespace zhurin_i_gauss_kernel_seq {

ZhurinIGaussKernelSEQ::ZhurinIGaussKernelSEQ(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = OutType{};
}

bool ZhurinIGaussKernelSEQ::ValidationImpl() {
  const auto &in = GetInput();
  int w = std::get<0>(in);
  int h = std::get<1>(in);
  int parts = std::get<2>(in);
  const auto &img = std::get<3>(in);

  if (w <= 0 || h <= 0 || parts <= 0 || parts > w) {
    return false;
  }
  if (static_cast<int>(img.size()) != h) {
    return false;
  }
  for (int i = 0; i < h; ++i) {
    if (static_cast<int>(img[i].size()) != w) {
      return false;
    }
  }
  return true;
}

bool ZhurinIGaussKernelSEQ::PreProcessingImpl() {
  const auto &in = GetInput();
  width = std::get<0>(in);
  height = std::get<1>(in);
  numParts = std::get<2>(in);
  image = std::get<3>(in);

  result.assign(height, std::vector<int>(width, 0));

  std::vector<std::vector<int>> padded(height + 2, std::vector<int>(width + 2, 0));
  for (int i = 0; i < height; ++i) {
    std::copy(image[i].begin(), image[i].end(), padded[i + 1].begin() + 1);
  }

  int baseWidth = width / numParts;
  int remainder = width % numParts;
  int xStart = 0;

  auto convolveAt = [&](int row, int col) -> int {
    int sum = 0;
    for (int ki = 0; ki < 3; ++ki) {
      for (int kj = 0; kj < 3; ++kj) {
        sum += padded[row - 1 + ki][col - 1 + kj] * kernel_[ki][kj];
      }
    }
    return sum >> shift_;
  };

  for (int part = 0; part < numParts; ++part) {
    int partWidth = baseWidth + (part < remainder ? 1 : 0);
    int xEnd = xStart + partWidth;

    for (int i = 1; i <= height; ++i) {
      for (int j = xStart + 1; j <= xEnd; ++j) {
        result[i - 1][j - 1] = convolveAt(i, j);
      }
    }

    xStart = xEnd;
  }

  GetOutput() = std::move(result);
  output_written = true;
  return true;
}

bool ZhurinIGaussKernelSEQ::RunImpl() {
  return true;
}

bool ZhurinIGaussKernelSEQ::PostProcessingImpl() {
  return output_written;
}

}  // namespace zhurin_i_gauss_kernel_seq
