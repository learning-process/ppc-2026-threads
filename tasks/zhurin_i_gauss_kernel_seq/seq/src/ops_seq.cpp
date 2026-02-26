#include "zhurin_i_gauss_kernel_seq/seq/include/ops_seq.hpp"

#include <algorithm>
#include <vector>

namespace zhurin_i_gauss_kernel_seq {

constexpr int ZhurinIGaussKernelSEQ::kernel[3][3];
constexpr int ZhurinIGaussKernelSEQ::slip;

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

  if (w <= 0 || h <= 0 || parts <= 0 || parts > w) return false;
  if (static_cast<int>(img.size()) != h) return false;
  for (int i = 0; i < h; ++i) {
    if (static_cast<int>(img[i].size()) != w) return false;
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
  return true;
}

bool ZhurinIGaussKernelSEQ::RunImpl() {
  if (width == 0 || height == 0 || numParts == 0) return false;

  std::vector<std::vector<int>> tmp(height + 2, std::vector<int>(width + 2, 0));

  for (int i = 0; i < height; ++i) {
    std::copy(image[i].begin(), image[i].end(), tmp[i + 1].begin() + 1);
  }

  int baseWidth = width / numParts;
  int remainder = width % numParts;

  int xStart = 0;
  for (int part = 0; part < numParts; ++part) {
    int partWidth = baseWidth + (part < remainder ? 1 : 0);
    int xEnd = xStart + partWidth;

    for (int i = 1; i <= height; ++i) {
      for (int j = xStart + 1; j <= xEnd; ++j) {
        int sum = 0;
        sum += tmp[i - 1][j - 1] * kernel[0][0];
        sum += tmp[i - 1][j]     * kernel[0][1];
        sum += tmp[i - 1][j + 1] * kernel[0][2];
        sum += tmp[i][j - 1]     * kernel[1][0];
        sum += tmp[i][j]         * kernel[1][1];
        sum += tmp[i][j + 1]     * kernel[1][2];
        sum += tmp[i + 1][j - 1] * kernel[2][0];
        sum += tmp[i + 1][j]     * kernel[2][1];
        sum += tmp[i + 1][j + 1] * kernel[2][2];

        result[i - 1][j - 1] = sum >> slip;  
      }
    }

    xStart = xEnd;
  }

  return true;
}

bool ZhurinIGaussKernelSEQ::PostProcessingImpl() {
  GetOutput() = std::move(result);
  return true;
}

}  // namespace zhurin_i_gauss_kernel_seq