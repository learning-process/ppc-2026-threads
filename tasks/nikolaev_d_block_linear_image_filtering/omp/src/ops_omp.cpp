#include "nikolaev_d_block_linear_image_filtering/omp/include/ops_omp.hpp"

#include <omp.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "nikolaev_d_block_linear_image_filtering/common/include/common.hpp"

namespace nikolaev_d_block_linear_image_filtering {

NikolaevDBlockLinearImageFilteringOMP::NikolaevDBlockLinearImageFilteringOMP(const InType &in) {
  SetTypeOfTask(GetStaticTypeOfTask());
  GetInput() = in;
  GetOutput() = std::vector<uint8_t>();
}

bool NikolaevDBlockLinearImageFilteringOMP::ValidationImpl() {
  const auto img_width = get<0>(GetInput());
  const auto img_height = get<1>(GetInput());
  const auto &pixel_data = get<2>(GetInput());

  return static_cast<std::size_t>(img_width) * static_cast<std::size_t>(img_height) * 3 == pixel_data.size();
}

bool NikolaevDBlockLinearImageFilteringOMP::PreProcessingImpl() {
  return true;
}

std::uint8_t NikolaevDBlockLinearImageFilteringOMP::GetPixel(const std::vector<uint8_t> &data, int w, int h, int nx,
                                                             int ny, int ch) {
  int ix = std::clamp(nx, 0, w - 1);
  int iy = std::clamp(ny, 0, h - 1);
  return data[((iy * w + ix) * 3) + ch];
}

bool NikolaevDBlockLinearImageFilteringOMP::RunImpl() {
  const int width = std::get<0>(GetInput());
  const int height = std::get<1>(GetInput());
  const auto &src = std::get<2>(GetInput());

  auto &dst = GetOutput();
  dst.assign(src.size(), 0);

  const std::array<std::array<int, 3>, 3> kernel = {{{1, 2, 1}, {2, 4, 2}, {1, 2, 1}}};
  const int sum = 16;

#pragma omp parallel for collapse(2) schedule(static) default(none) shared(width, height, src, dst, kernel, sum)
  for (int ny = 0; ny < height; ++ny) {
    for (int nx = 0; nx < width; ++nx) {
      for (int ch = 0; ch < 3; ++ch) {
        int acc = 0;
        for (int ky = -1; ky <= 1; ++ky) {
          for (int kx = -1; kx <= 1; ++kx) {
            acc += GetPixel(src, width, height, nx + kx, ny + ky, ch) * kernel.at(ky + 1).at(kx + 1);
          }
        }

        int res = (acc + 8) / sum;
        dst[((ny * width + nx) * 3) + ch] = static_cast<uint8_t>(std::clamp(res, 0, 255));
      }
    }
  }
  return true;
}

bool NikolaevDBlockLinearImageFilteringOMP::PostProcessingImpl() {
  return true;
}

}  // namespace nikolaev_d_block_linear_image_filtering
